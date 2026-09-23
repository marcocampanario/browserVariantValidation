#!/usr/bin/env Rscript
# Validate HGVS annotations and compare their genomic mapping with a browser TSV.
# Adapted from the validation logic in marcocampanario/hgvs2vcf (2026-09-23).
# Standalone
# Dependency for online requests: install.packages("httr2")
#
# R:
#   source("validar_browser.R")
#   triagem <- validar_browser("browser.tsv", assembly = "GRCh38", consultar = FALSE)
#   resultado <- validar_browser("browser.tsv", assembly = "GRCh38")
#   salvar_browser_tsv(resultado, "browser_validado.tsv")
#
# Terminal (the assembly is REQUIRED):
#   Rscript validar_browser.R browser.tsv browser_validado.tsv GRCh38
#   Rscript validar_browser.R browser.tsv browser_triagem.tsv GRCh38 --triagem
#
# Required columns: CHROM POS REF ALT TRANSCRIPT HGVS_C.
# All original columns, values, row order and duplicates are retained. Only
# additional columns are written. Input TSV fields are read as character strings
# with check.names = FALSE, preserving names such as FREQ_CENTRO-OESTE.
# POS must be 1-based; use one explicit REF/ALT allele and one annotation per row.
# TRANSCRIPT + HGVS_C form the query; a complete HGVS_C is also accepted.
# Versioned RefSeq references and human Ensembl transcripts (ENST...) are accepted.
# Missing versions, multiple annotations and conflicting transcript IDs are
# flagged, never guessed. No transcript is inferred from GENE or from other rows.
#
# Interpretation:
#   VV_STATUS: HGVS/API result (OK, REVIEW, API_ERROR, READY in screening, etc.).
#   VCF_INPUT_STATUS: whether the input VCF fields can be compared.
#   CHECK_CHROM/POS/REF/ALT: TRUE/FALSE; NA means not evaluated.
#   VCF_COMPARISON: MATCH, SNV_MISMATCH, REPRESENTATION_DIFF_REVIEW,
#                   INPUT_VCF_REVIEW, or NOT_EVALUATED.
#   VALIDATION_STATUS: PASS only when VV_STATUS == OK AND VCF_COMPARISON == MATCH;
#                      REVIEW for a mapped result needing inspection;
#                      NOT_EVALUATED when no usable mapping was obtained.
#
# Limits:
# - This checks consistency between DNA HGVS and genomic fields. It does not
#   independently check REF against a local FASTA, normalize/left-align the input,
#   perform liftover, or validate frequencies, genotypes, IMPACT or EFFECT.
# - Different INDEL/MNV representations require normalization against the SAME
#   assembly FASTA before biological inequivalence can be concluded.
# - HGVS_P is retained; VV_HGVS_P is the predicted result, not a validation of
#   the supplied protein annotation. Gene symbols are returned for inspection.
# - Missing HGVS has no coordinate-only fallback: it remains NOT_EVALUATED.
# - API warnings/corrections never silently become PASS.
# - Each distinct HGVS is requested once per run. Public-service throughput is
#   unsuitable for unrestricted whole-exome/genome tables; a local/authorized
#   endpoint can be provided via base_url in R or VV_BASE_URL in the terminal.
# - Only the assembled DNA HGVS and assembly are sent to VariantValidator.
#   If required by your endpoint, set an authorized bearer token via VV_TOKEN.
# - API metadata and retrieval time are retained in attributes(result); saveRDS
#   can preserve these in addition to the TSV. A TSV contains VV_ASSEMBLY.
#
# API documentation:
# https://openvar.github.io/variantValidator/rest-vv/rest_VariantValidator.html

.hb_text <- function(x) {
  if (is.null(x) || !is.atomic(x) || length(x) != 1L || is.na(x)) return(NA_character_)
  as.character(x)
}

.hb_note <- function(...) {
  x <- as.character(unlist(list(...), use.names = FALSE))
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(NA_character_)
  gsub("[\t\r\n]+", " ", paste(x, collapse = " | "))
}

.hb_missing <- function(x) is.na(x) | trimws(x) %in% c("", ".", "NA")

.hb_chrom <- function(x) {
  x <- toupper(sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE))
  x[x %in% c("M", "MT")] <- "MT"
  x
}

.hb_query <- function(transcript, hgvs) {
  fail <- function(status, note) list(query = NA_character_, status = status, note = note)
  transcript <- trimws(as.character(transcript))
  hgvs <- trimws(as.character(hgvs))
  if (.hb_missing(hgvs)) return(fail("MISSING_HGVS", "HGVS_C ausente; nao foi inferido a partir do VCF."))
  if (grepl("[;,|]", hgvs) || (!.hb_missing(transcript) && grepl("[;,|]", transcript))) {
    return(fail("MULTIPLE_ANNOTATIONS", "Use uma anotacao/transcrito por linha; nenhuma lista foi truncada."))
  }
  if (grepl("[\\[\\]]", hgvs, perl = TRUE)) {
    return(fail("COMPLEX_HGVS_REVIEW", "Alelos complexos exigem tratamento especifico."))
  }
  # Accept a full HGVS string, or combine a bare c./n. description with TRANSCRIPT.
  if (grepl(":", hgvs, fixed = TRUE)) {
    parts <- strsplit(hgvs, ":", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) return(fail("INVALID_HGVS_FORMAT", "Esperado accession.version:c./n./g./m...."))
    accession <- trimws(parts[[1L]])
    change <- trimws(parts[[2L]])
    if (grepl("^(NM_|NR_|ENST)", accession) && !.hb_missing(transcript) &&
        !identical(transcript, accession)) {
      return(fail("TRANSCRIPT_CONFLICT", "TRANSCRIPT difere da referencia em HGVS_C; confira a versao."))
    }
  } else {
    if (.hb_missing(transcript)) return(fail("MISSING_TRANSCRIPT", "Informe TRANSCRIPT com accession.version."))
    accession <- transcript
    change <- hgvs
  }
  if (!grepl("^(?:(?:NM|NR|NC|NG|NT|NW)_[0-9]+|ENST[0-9]+)\\.[0-9]+$", accession, perl = TRUE)) {
    status <- if (grepl("^(?:(?:NM|NR|NC|NG|NT|NW)_[0-9]+|ENST[0-9]+)$", accession, perl = TRUE))
      "MISSING_REFERENCE_VERSION" else "UNSUPPORTED_REFERENCE"
    return(fail(status, "Informe uma referencia RefSeq ou ENST com versao; nenhuma versao foi presumida."))
  }
  if (!grepl("^[cngm]\\.\\S+$", change, perl = TRUE)) {
    return(fail("INVALID_HGVS_FORMAT", "Esperada uma descricao de DNA c., n., g. ou m. sem espacos internos."))
  }
  query <- paste0(accession, ":", change)
  # Large genomic deletions need SV-aware output, as in the original program.
  deletion <- regmatches(change, regexec("^g\\.([0-9]+)(?:_([0-9]+))?del(?:[ACGT]*)$", change, perl = TRUE))[[1L]]
  if (length(deletion) == 3L) {
    start <- as.numeric(deletion[[2L]])
    end <- if (nzchar(deletion[[3L]])) as.numeric(deletion[[3L]]) else start
    if (start < 1 || end < start) return(fail("INVALID_INTERVAL", "Intervalo de delecao invalido."))
    if (end - start + 1 >= 50) return(fail("SV_REVIEW", "Delecao genomica >=50 bp; requer representacao SV/END."))
  }
  list(query = query, status = "READY", note = NA_character_)
}

.hb_request <- function(query, assembly, base_url, token) {
  accession <- sub(":.*$", "", query)
  endpoint <- if (startsWith(accession, "ENST")) "variantvalidator_ensembl" else "variantvalidator"
  select <- if (grepl("^(NM_|NR_|ENST)", accession)) accession else "mane_select"
  components <- c("VariantValidator", endpoint, assembly, query, select)
  encoded <- vapply(components, utils::URLencode, character(1L), reserved = TRUE)
  url <- paste0(sub("/+$", "", base_url), "/", paste(encoded, collapse = "/"))
  req <- httr2::request(url)
  req <- httr2::req_url_query(req, `content-type` = "application/json")
  req <- httr2::req_headers(req, Accept = "application/json")
  req <- httr2::req_user_agent(req, "hgvs2vcf-browser-R/1.0")
  req <- httr2::req_timeout(req, 90)
  req <- httr2::req_retry(req, max_tries = 3L)
  if (nzchar(token)) req <- httr2::req_auth_bearer_token(req, token)
  httr2::resp_body_json(httr2::req_perform(req), simplifyVector = FALSE)
}

.hb_parse <- function(data, query, assembly) {
  if (!is.list(data) || is.null(names(data))) stop("Resposta da API nao e um objeto JSON nomeado.")
  nodes <- Filter(is.list, data[setdiff(names(data), c("metadata", "flag"))])
  warnings <- .hb_note(lapply(nodes, function(node) node$validation_warnings))
  fail <- function(status, note) list(VV_STATUS = status, VV_WARNINGS = .hb_note(note, warnings))
  candidates <- Filter(function(node) {
    identical(.hb_text(node$submitted_variant), query) &&
      is.list(node$primary_assembly_loci[[tolower(assembly)]])
  }, nodes)
  if (!length(candidates)) return(fail("NO_TARGET_MAPPING", "Sem mapeamento explicito para a montagem alvo."))
  accession <- sub(":.*$", "", query)
  is_transcript <- grepl("^(NM_|NR_|ENST)", accession)
  if (is_transcript) {
    candidates <- Filter(function(node) {
      returned <- .hb_text(node$hgvs_transcript_variant)
      !is.na(returned) && identical(sub(":.*$", "", returned), accession)
    }, candidates)
    if (!length(candidates)) return(fail("TRANSCRIPT_VERSION_MISMATCH", "A API nao preservou o transcrito e sua versao."))
  }
  signatures <- vapply(candidates, function(node) {
    v <- node$primary_assembly_loci[[tolower(assembly)]]$vcf
    paste(.hb_text(v$chr), .hb_text(v$pos), .hb_text(v$ref), .hb_text(v$alt), sep = ":")
  }, character(1L))
  if (length(unique(signatures)) != 1L) return(fail("AMBIGUOUS_MAPPING", "Mais de uma representacao genomica foi retornada."))
  node <- candidates[[1L]]
  locus <- node$primary_assembly_loci[[tolower(assembly)]]
  vcf <- locus$vcf
  chr <- .hb_chrom(.hb_text(vcf$chr))
  pos_text <- .hb_text(vcf$pos)
  pos <- suppressWarnings(as.numeric(pos_text))
  ref <- toupper(.hb_text(vcf$ref))
  alt <- toupper(.hb_text(vcf$alt))
  if (anyNA(c(chr, pos, ref, alt)) || !is.finite(pos) || pos < 1 || pos != floor(pos) ||
      !grepl("^[ACGTN]+$", ref) || !grepl("^[ACGTN]+$", alt) || identical(ref, alt)) {
    return(fail("UNSUPPORTED_VCF_ALLELES", "Resposta sem variante VCF explicita e comparavel."))
  }
  if (!chr %in% c(as.character(1:22), "X", "Y", "MT")) {
    return(fail("NON_PRIMARY_CONTIG", "Mapeamento fora dos cromossomos primarios esperados."))
  }
  hgvs_c <- .hb_text(node$hgvs_transcript_variant)
  hgvs_g <- .hb_text(locus$hgvs_genomic_description)
  returned <- if (is_transcript) hgvs_c else hgvs_g
  if (!identical(returned, query)) warnings <- .hb_note(warnings, "HGVS retornado difere da consulta; revisar correcao/normalizacao/remapeamento.")
  if (grepl("N", paste0(ref, alt))) warnings <- .hb_note(warnings, "Alelo retornado contem base N.")
  list(VV_STATUS = if (is.na(warnings)) "OK" else "REVIEW", VV_WARNINGS = warnings,
       VV_CHROM = chr, VV_POS = pos, VV_REF = ref, VV_ALT = alt,
       VV_GENE = .hb_text(node$gene_symbol), VV_HGVS_C = hgvs_c, VV_HGVS_G = hgvs_g,
       VV_HGVS_P = .hb_text(node$hgvs_predicted_protein_consequence$tlr))
}

.hb_read <- function(path) {
  if (length(path) != 1L || !file.exists(path) || dir.exists(path)) stop("TSV de entrada nao encontrado.", call. = FALSE)
  fields <- utils::count.fields(path, sep = "\t", quote = "", comment.char = "", blank.lines.skip = TRUE)
  if (!length(fields) || anyNA(fields) || any(fields != fields[[1L]])) {
    stop("TSV irregular: todas as linhas devem ter o mesmo numero de campos do cabecalho.", call. = FALSE)
  }
  utils::read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                    colClasses = "character", na.strings = NULL, check.names = FALSE,
                    stringsAsFactors = FALSE, row.names = NULL, fill = FALSE,
                    fileEncoding = "UTF-8-BOM")
}

validar_browser <- function(tabela, assembly, consultar = TRUE,
                            base_url = Sys.getenv("VV_BASE_URL", "https://rest.variantvalidator.org"),
                            token = Sys.getenv("VV_TOKEN", ""), pausa = 0.4) {
  if (missing(assembly) || length(assembly) != 1L || is.na(assembly) ||
      !assembly %in% c("GRCh38", "GRCh37")) {
    stop("Informe assembly = 'GRCh38' ou 'GRCh37', conforme a montagem da tabela.", call. = FALSE)
  }
  if (!is.logical(consultar) || length(consultar) != 1L || is.na(consultar)) stop("consultar deve ser TRUE ou FALSE.", call. = FALSE)
  if (!is.numeric(pausa) || length(pausa) != 1L || !is.finite(pausa) || pausa < 0) stop("pausa deve ser >=0.", call. = FALSE)
  if (is.character(tabela) && length(tabela) == 1L) tabela <- .hb_read(tabela)
  if (!is.data.frame(tabela)) stop("tabela deve ser um data.frame ou caminho de TSV.", call. = FALSE)
  tabela <- as.data.frame(tabela, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(names(tabela))) stop("Ha nomes de colunas duplicados na entrada.", call. = FALSE)
  required <- c("CHROM", "POS", "REF", "ALT", "TRANSCRIPT", "HGVS_C")
  absent <- setdiff(required, names(tabela))
  if (length(absent)) stop("Colunas ausentes: ", paste(absent, collapse = ", "), call. = FALSE)
  if (any(!vapply(tabela, function(x) is.atomic(x) && is.null(dim(x)), logical(1L)))) {
    stop("Use colunas simples, sem listas ou matrizes.", call. = FALSE)
  }
  n <- nrow(tabela)
  extra <- data.frame(
    VV_ROW_ID = seq_len(n), VV_ASSEMBLY = rep(assembly, n),
    VV_QUERY = rep(NA_character_, n), VV_STATUS = rep(NA_character_, n),
    VV_CHROM = rep(NA_character_, n), VV_POS = rep(NA_real_, n),
    VV_REF = rep(NA_character_, n), VV_ALT = rep(NA_character_, n),
    VV_GENE = rep(NA_character_, n), VV_HGVS_C = rep(NA_character_, n),
    VV_HGVS_G = rep(NA_character_, n), VV_HGVS_P = rep(NA_character_, n),
    VV_WARNINGS = rep(NA_character_, n), stringsAsFactors = FALSE
  )
  added <- c(names(extra), "VCF_INPUT_STATUS", "CHECK_CHROM", "CHECK_POS", "CHECK_REF", "CHECK_ALT", "VCF_COMPARISON", "VALIDATION_STATUS")
  collision <- intersect(names(tabela), added)
  if (length(collision)) stop("A entrada ja contem colunas de validacao: ", paste(collision, collapse = ", "), call. = FALSE)
  for (i in seq_len(n)) {
    info <- .hb_query(tabela$TRANSCRIPT[[i]], tabela$HGVS_C[[i]])
    extra$VV_QUERY[[i]] <- info$query
    extra$VV_STATUS[[i]] <- info$status
    extra$VV_WARNINGS[[i]] <- info$note
  }
  ready <- which(extra$VV_STATUS == "READY")
  queries <- unique(extra$VV_QUERY[ready])
  metadata <- list()
  if (consultar && length(queries)) {
    if (!requireNamespace("httr2", quietly = TRUE)) stop("Instale httr2: install.packages('httr2').", call. = FALSE)
    # Query each distinct description once, then map back by row index (no join).
    results <- vector("list", length(queries))
    for (j in seq_along(queries)) {
      if (j == 1L || j %% 25L == 0L || j == length(queries)) message("Consultando HGVS unico ", j, "/", length(queries), ".")
      if (pausa > 0) Sys.sleep(pausa)
      data <- tryCatch(.hb_request(queries[[j]], assembly, base_url, token), error = identity)
      if (inherits(data, c("httr2_http_401", "httr2_http_403"))) {
        stop("A API recusou acesso. Verifique as credenciais autorizadas (VV_TOKEN) e o endpoint.", call. = FALSE)
      }
      if (inherits(data, "error")) {
        results[[j]] <- list(VV_STATUS = "API_ERROR", VV_WARNINGS = .hb_note(conditionMessage(data)))
      } else {
        if (is.list(data) && !is.null(data$metadata) && !any(vapply(metadata, identical, logical(1L), data$metadata))) metadata[[length(metadata) + 1L]] <- data$metadata
        results[[j]] <- tryCatch(.hb_parse(data, queries[[j]], assembly), error = function(e) {
          list(VV_STATUS = "RESPONSE_ERROR", VV_WARNINGS = .hb_note(conditionMessage(e)))
        })
      }
    }
    index <- match(extra$VV_QUERY[ready], queries)
    for (k in seq_along(ready)) {
      i <- ready[[k]]
      parsed <- results[[index[[k]]]]
      for (name in names(parsed)) extra[[name]][[i]] <- parsed[[name]]
    }
  }

  # Compare genomic representations in the same assembly and on the forward strand.
  # Chromosome aliases and letter case are harmonized only in comparison vectors.
  chr <- .hb_chrom(tabela$CHROM)
  pos_text <- trimws(as.character(tabela$POS))
  pos <- suppressWarnings(as.numeric(pos_text))
  ref <- toupper(trimws(as.character(tabela$REF)))
  alt <- toupper(trimws(as.character(tabela$ALT)))
  input_status <- rep("OK", n)
  input_status[!(chr %in% c(as.character(1:22), "X", "Y", "MT"))] <- "NON_PRIMARY_CONTIG"
  input_status[!grepl("^[ACGTN]+$", ref) | !grepl("^[ACGTN]+$", alt)] <- "UNSUPPORTED_ALLELES"
  input_status[grepl("N", ref) | grepl("N", alt)] <- "AMBIGUOUS_BASES"
  input_status[grepl(",", alt, fixed = TRUE)] <- "MULTIALLELIC"
  input_status[which(!is.na(ref) & !is.na(alt) & ref == alt)] <- "SAME_REF_ALT"
  input_status[!grepl("^[0-9]+$", pos_text) | is.na(pos) | !is.finite(pos) | pos < 1 | pos > .Machine$integer.max] <- "INVALID_POS"
  missing_fields <- Reduce(`|`, lapply(tabela[c("CHROM", "POS", "REF", "ALT")], function(x) .hb_missing(as.character(x))))
  input_status[missing_fields] <- "MISSING_FIELDS"
  extra$VCF_INPUT_STATUS <- input_status
  mapped <- extra$VV_STATUS %in% c("OK", "REVIEW")
  comparable <- mapped & input_status == "OK"
  extra$CHECK_CHROM <- extra$CHECK_POS <- extra$CHECK_REF <- extra$CHECK_ALT <- rep(NA, n)
  extra$CHECK_CHROM[comparable] <- chr[comparable] == extra$VV_CHROM[comparable]
  extra$CHECK_POS[comparable] <- pos[comparable] == extra$VV_POS[comparable]
  extra$CHECK_REF[comparable] <- ref[comparable] == extra$VV_REF[comparable]
  extra$CHECK_ALT[comparable] <- alt[comparable] == extra$VV_ALT[comparable]
  same <- comparable & Reduce(`&`, extra[c("CHECK_CHROM", "CHECK_POS", "CHECK_REF", "CHECK_ALT")])
  snvs <- nchar(ref) == 1L & nchar(alt) == 1L & nchar(extra$VV_REF) == 1L & nchar(extra$VV_ALT) == 1L
  extra$VCF_COMPARISON <- rep("NOT_EVALUATED", n)
  extra$VCF_COMPARISON[mapped & !comparable] <- "INPUT_VCF_REVIEW"
  extra$VCF_COMPARISON[comparable] <- "REPRESENTATION_DIFF_REVIEW"
  extra$VCF_COMPARISON[which(comparable & snvs)] <- "SNV_MISMATCH"
  extra$VCF_COMPARISON[which(same)] <- "MATCH"
  extra$VALIDATION_STATUS <- rep("NOT_EVALUATED", n)
  extra$VALIDATION_STATUS[mapped] <- "REVIEW"
  extra$VALIDATION_STATUS[extra$VV_STATUS == "OK" & extra$VCF_COMPARISON == "MATCH"] <- "PASS"
  result <- cbind(tabela, extra)
  attr(result, "vv_metadata") <- metadata
  attr(result, "retrieval_time_utc") <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  attr(result, "assembly") <- assembly
  result
}

salvar_browser_tsv <- function(tabela, arquivo) {
  if (file.exists(arquivo)) stop("A saida ja existe: ", arquivo, call. = FALSE)
  if (!dir.exists(dirname(arquivo))) stop("A pasta de saida nao existe.", call. = FALSE)
  # Unquoted TSV: reject embedded separators instead of silently damaging fields.
  if (any(grepl("[\t\r\n]", names(tabela))) || any(vapply(tabela, function(x) any(grepl("[\t\r\n]", as.character(x))), logical(1L)))) {
    stop("Ha campos com tabulacao/quebra de linha; revise antes de exportar TSV.", call. = FALSE)
  }
  utils::write.table(tabela, file = arquivo, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
  invisible(arquivo)
}

.hb_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste("Uso: Rscript validar_browser.R entrada.tsv saida.tsv GRCh38|GRCh37 [--triagem]",
                 "A montagem deve ser a mesma usada na tabela; POS e 1-based.", sep = "\n")
  if (any(args %in% c("-h", "--help"))) { cat(usage, "\n"); return(invisible(NULL)) }
  screening <- "--triagem" %in% args
  args <- args[args != "--triagem"]
  if (length(args) != 3L || any(startsWith(args, "--"))) stop(usage, call. = FALSE)
  if (file.exists(args[[2L]])) stop("A saida ja existe: ", args[[2L]], call. = FALSE)
  if (!dir.exists(dirname(args[[2L]]))) stop("A pasta de saida nao existe.", call. = FALSE)
  result <- validar_browser(args[[1L]], assembly = args[[3L]], consultar = !screening)
  salvar_browser_tsv(result, args[[2L]])
  print(table(result$VV_STATUS, useNA = "ifany"))
  print(table(result$VCF_COMPARISON, useNA = "ifany"))
  message("Tabela salva: ", args[[2L]])
  invisible(result)
}

if (sys.nframe() == 0L) {
  tryCatch(.hb_main(), error = function(e) {
    message("ERRO: ", conditionMessage(e))
    quit(save = "no", status = 1L)
  })
}

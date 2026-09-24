#!/usr/bin/env Rscript
# Validate a variant table against the NCBI dbSNP reference VCF.
# Standalone. Needs one of: Bioconductor package Rsamtools, or `tabix` (htslib) in PATH.
#   install.packages("BiocManager"); BiocManager::install("Rsamtools")
#
# Terminal use (assembly and --dbsnp are required):
#   Rscript validate_dbsnp.R [input.tsv] [output.tsv] [assembly] --dbsnp [path/to/dbsnp.gz]
#   Writes: full comparison table (table.dbsnp.tsv) and problems table (table_dbsnp_problemas.tsv)
#   Options: --scan-rsids          search the whole VCF for rsIDs not found at POS (slow)
#            --janela [N]            bp tolerance to find the same rsID nearby (default 10)
#            --vocabulario [ARQ]     extra TSV (termo<TAB>flags) added to the built-in vocabulary
#            --col [TYPE=COLUMN]    force a column, e.g. --col gene=GENE_SYMBOL (repeatable)
#            --problemas [ARQ]      name of the problems file (default: <saida>_problemas.tsv)
#
# Columns (found by name, case-insensitive; override with --col / colunas = list(...)):
#   required: CHROM, POS
#   optional: REF, ALT, RSID (or ID/SNP/...), GENE, REGION, EFFECT, IMPACT
# Checks whose column is absent are simply not evaluated (NA).
# All original columns, values, row order and duplicates are retained. Only
# additional columns are written. Input fields are read as character strings.
#
# What is compared, only against dbSNP (no FASTA, no format/syntax report):
#   CHECK_RSID     rsID is the dbSNP record at CHROM:POS (+-janela bp for indels)
#   CHECK_ALLELES  REF/ALT match the dbSNP record (minimal representation; "-" = empty allele)
#   CHECK_GENE     every gene in the row is in dbSNP GENEINFO/PSEUDOGENEINFO
#   CHECK_REGION / CHECK_EFFECT / CHECK_IMPACT
#                  every term is supported by the dbSNP functional flags
#                  (NSM NSN NSF SYN ASS DSS U5 U3 INT R5 R3; no flag + no gene = intergenic)
#
# Interpretation:
#   DBSNP_STATUS: CONFIRMED, NEIGHBOR_POSITION, RSID_MISMATCH, NOT_AT_POSITION,
#                 RSID_ELSEWHERE, RSID_NOT_IN_DBSNP (these two need --scan-rsids),
#                 RSID_AVAILABLE (row has no rsID, dbSNP has one), NOVEL (nothing in
#                 dbSNP), CONTIG_NOT_IN_DBSNP, NOT_EVALUATED (unreadable CHROM/POS).
#   CHECK_*: TRUE/FALSE; NA means not evaluated.
#   DBSNP_VALIDATION: FAIL if any CHECK_* is FALSE or the rsID does not match;
#                     REVIEW if something needs inspection (see DBSNP_NOTES);
#                     PASS if everything evaluated agrees; NOT_EVALUATED otherwise.
#
# Problems file (one row per problem; a row with two problems appears twice):
#   ROW (data row of the input, 1 = first after the header), CHROM, POS, ID, REF, ALT,
#   GENE, DBSNP_VALIDATION, LEVEL (FAIL/REVIEW), TYPE, MESSAGE, DBSNP_RSID, DBSNP_POS.
#   TYPE: RSID_MISMATCH, RSID_NOT_AT_POSITION, RSID_ELSEWHERE, RSID_NOT_IN_DBSNP,
#         NEIGHBOR_POSITION, RSID_AVAILABLE, ALLELES_MISMATCH, ALLELES_REVERSE_STRAND,
#         ALT_NOT_IN_DBSNP, GENE_MISMATCH, GENE_EMPTY, GENE_ROW_MISSING,
#         {REGION,EFFECT,IMPACT}_MISMATCH / _UNKNOWN_TERM / _NOT_VERIFIABLE,
#         EFFECT_IMPACT_INCONSISTENT.
#
# Conventions:
# - A variant in two genes has one row per gene: rows linked to the same dbSNP
#   record are checked together, and a dbSNP gene with no row is noted (REVIEW).
# - IMPACT (HIGH/MODERATE/LOW/MODIFIER) is also checked against EFFECT in the same row.
#
# Limits (from dbSNP itself):
# - Functional flags are per rsID, not per allele, transcript or gene.
# - No flags exist for inframe indels, start/stop lost, splice region or non-coding
#   exons: those terms are recognised but not checked.
# - dbSNP links genes within 2 kb upstream / 500 bp downstream (SnpEff uses 5 kb).
# - Merged/retired rsIDs are not in the VCF.
#
# Constants for compability
.dv_FLAGS <- c(NSF = "frameshift_variant", NSM = "missense_variant", NSN = "stop_gained",
               SYN = "synonymous_variant", ASS = "splice_acceptor_variant",
               DSS = "splice_donor_variant", U5 = "5_prime_UTR_variant",
               U3 = "3_prime_UTR_variant", INT = "intron_variant",
               R5 = "upstream_gene_variant", R3 = "downstream_gene_variant")
.dv_REGION <- c(NSF = "coding", NSM = "coding", NSN = "coding", SYN = "coding",
                ASS = "splice_site", DSS = "splice_site", U5 = "5'UTR", U3 = "3'UTR",
                INT = "intron", R5 = "upstream", R3 = "downstream")
.dv_CLASS <- c(NSF = "HIGH", NSN = "HIGH", ASS = "HIGH", DSS = "HIGH", NSM = "MODERATE",
               SYN = "LOW", INT = "MODIFIER", U3 = "MODIFIER", U5 = "MODIFIER",
               R3 = "MODIFIER", R5 = "MODIFIER", INTERGENIC = "MODIFIER")
.dv_CLASS_TERMS <- c(high = "HIGH", alto = "HIGH", moderate = "MODERATE",
                     moderado = "MODERATE", low = "LOW", baixo = "LOW",
                     modifier = "MODIFIER", modificador = "MODIFIER")
.dv_UNCHECKABLE <- "NAO_VERIFICAVEL"
.dv_EXPECTED_INFO <- c("GENEINFO", "VC", names(.dv_FLAGS))

.dv_ALIASES <- list(
  chrom  = c("chrom", "chr", "chromosome", "chromosome_name", "chr_name"),
  pos    = c("pos", "position", "bp", "base_pair_location", "start"),
  ref    = c("ref", "reference", "ref_allele", "reference_allele"),
  alt    = c("alt", "alternate", "alt_allele", "alternate_allele"),
  rsid   = c("rsid", "rs_id", "dbsnp", "dbsnp_id", "rs", "snp", "snp_id", "variant_id", "id"),
  gene   = c("gene", "genes", "gene_symbol", "symbol", "hgnc_symbol", "gene_name",
             "gene.refgene"),
  region = c("region", "regiao", "gene_region", "genomic_region", "func.refgene", "func"),
  effect = c("effect", "efeito", "consequence", "consequencia", "annotation",
             "variant_effect", "most_severe_consequence", "exonicfunc.refgene"),
  impact = c("impact", "impacto", "putative_impact", "annotation_impact", "impact_class")
)

# Term -> dbSNP flags that confirm it. INTERGENIC = no flag and no gene.
# NAO_VERIFICAVEL = recognised term without a dbSNP flag (not checked, no note).
.dv_VOCAB_TEXT <- c(
  "coding\tNSF,NSM,NSN,SYN",
  "coding_sequence_variant\tNSF,NSM,NSN,SYN",
  "cds\tNSF,NSM,NSN,SYN",
  "exonic\tNSF,NSM,NSN,SYN",
  "exon\tNSF,NSM,NSN,SYN",
  "exonico\tNSF,NSM,NSN,SYN",
  "exonica\tNSF,NSM,NSN,SYN",
  "codificante\tNSF,NSM,NSN,SYN",
  "utr\tU3,U5",
  "utr_variant\tU3,U5",
  "utr5\tU5",
  "5_utr\tU5",
  "5utr\tU5",
  "utr_5\tU5",
  "5_prime_utr\tU5",
  "5_prime_utr_variant\tU5",
  "utr3\tU3",
  "3_utr\tU3",
  "3utr\tU3",
  "utr_3\tU3",
  "3_prime_utr\tU3",
  "3_prime_utr_variant\tU3",
  "intron\tINT",
  "intronic\tINT",
  "intron_variant\tINT",
  "intronico\tINT",
  "intronica\tINT",
  "splicing\tASS,DSS",
  "splice\tASS,DSS",
  "splice_site\tASS,DSS",
  "splice_site_variant\tASS,DSS",
  "upstream\tR5",
  "upstream_gene_variant\tR5",
  "upstream_variant\tR5",
  "2kb_upstream_variant\tR5",
  "downstream\tR3",
  "downstream_gene_variant\tR3",
  "downstream_variant\tR3",
  "500b_downstream_variant\tR3",
  "intergenic\tINTERGENIC",
  "intergenic_variant\tINTERGENIC",
  "intergenic_region\tINTERGENIC",
  "intergenico\tINTERGENIC",
  "intergenica\tINTERGENIC",
  "missense\tNSM",
  "missense_variant\tNSM",
  "missense_mutation\tNSM",
  "nonsynonymous\tNSM",
  "non_synonymous\tNSM",
  "nonsynonymous_snv\tNSM",
  "non_synonymous_coding\tNSM",
  "nsm\tNSM",
  "stop_gained\tNSN",
  "stopgain\tNSN",
  "stop_gain\tNSN",
  "nonsense\tNSN",
  "nonsense_mutation\tNSN",
  "nsn\tNSN",
  "frameshift\tNSF",
  "frameshift_variant\tNSF",
  "frameshift_deletion\tNSF",
  "frameshift_insertion\tNSF",
  "frameshift_substitution\tNSF",
  "frame_shift\tNSF",
  "nsf\tNSF",
  "synonymous\tSYN",
  "synonymous_variant\tSYN",
  "synonymous_snv\tSYN",
  "silent\tSYN",
  "sinonima\tSYN",
  "sinonimo\tSYN",
  "syn\tSYN",
  "splice_acceptor_variant\tASS",
  "splice_acceptor\tASS",
  "ass\tASS",
  "splice_donor_variant\tDSS",
  "splice_donor\tDSS",
  "dss\tDSS",
  "int\tINT",
  "u3\tU3",
  "u5\tU5",
  "r3\tR3",
  "r5\tR5",
  "high\tNSF,NSN,ASS,DSS",
  "alto\tNSF,NSN,ASS,DSS",
  "moderate\tNSM",
  "moderado\tNSM",
  "low\tSYN",
  "baixo\tSYN",
  "modifier\tINT,U3,U5,R3,R5,INTERGENIC",
  "modificador\tINT,U3,U5,R3,R5,INTERGENIC",
  "intragenic_variant\tNAO_VERIFICAVEL",
  "5_prime_utr_premature_start_codon_gain_variant\tU5",
  "utr_5_prime\tU5",
  "utr_3_prime\tU3",
  "splice_region_variant\tNAO_VERIFICAVEL",
  "splice_donor_5th_base_variant\tNAO_VERIFICAVEL",
  "splice_donor_region_variant\tNAO_VERIFICAVEL",
  "splice_polypyrimidine_tract_variant\tNAO_VERIFICAVEL",
  "non_coding_transcript_exon_variant\tNAO_VERIFICAVEL",
  "non_coding_exon_variant\tNAO_VERIFICAVEL",
  "non_coding_transcript_variant\tNAO_VERIFICAVEL",
  "nmd_transcript_variant\tNAO_VERIFICAVEL",
  "mature_mirna_variant\tNAO_VERIFICAVEL",
  "inframe_insertion\tNAO_VERIFICAVEL",
  "inframe_deletion\tNAO_VERIFICAVEL",
  "conservative_inframe_insertion\tNAO_VERIFICAVEL",
  "conservative_inframe_deletion\tNAO_VERIFICAVEL",
  "disruptive_inframe_insertion\tNAO_VERIFICAVEL",
  "disruptive_inframe_deletion\tNAO_VERIFICAVEL",
  "nonframeshift_deletion\tNAO_VERIFICAVEL",
  "nonframeshift_insertion\tNAO_VERIFICAVEL",
  "nonframeshift_substitution\tNAO_VERIFICAVEL",
  "start_lost\tNAO_VERIFICAVEL",
  "startloss\tNAO_VERIFICAVEL",
  "initiator_codon_variant\tNAO_VERIFICAVEL",
  "stop_lost\tNAO_VERIFICAVEL",
  "stoploss\tNAO_VERIFICAVEL",
  "stop_retained_variant\tNAO_VERIFICAVEL",
  "start_retained_variant\tNAO_VERIFICAVEL",
  "protein_altering_variant\tNAO_VERIFICAVEL",
  "incomplete_terminal_codon_variant\tNAO_VERIFICAVEL",
  "coding_transcript_variant\tNAO_VERIFICAVEL",
  "exon_loss_variant\tNAO_VERIFICAVEL",
  "gene_fusion\tNAO_VERIFICAVEL",
  "bidirectional_gene_fusion\tNAO_VERIFICAVEL",
  "feature_ablation\tNAO_VERIFICAVEL",
  "transcript_ablation\tNAO_VERIFICAVEL",
  "transcript_amplification\tNAO_VERIFICAVEL",
  "sequence_feature\tNAO_VERIFICAVEL",
  "structural_interaction_variant\tNAO_VERIFICAVEL",
  "protein_protein_contact\tNAO_VERIFICAVEL",
  "tf_binding_site_variant\tNAO_VERIFICAVEL",
  "tfbs_ablation\tNAO_VERIFICAVEL",
  "regulatory_region_variant\tNAO_VERIFICAVEL",
  "regulatory_region_ablation\tNAO_VERIFICAVEL",
  "ncrna_exonic\tNAO_VERIFICAVEL",
  "ncrna_intronic\tINT",
  "ncrna_splicing\tASS,DSS",
  "unknown\tNAO_VERIFICAVEL"
)

# ---- small helpers --------------------------------------------------------------

.dv_note <- function(...) {
  x <- as.character(unlist(list(...), use.names = FALSE))
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(NA_character_)
  gsub("[\t\r\n]+", " ", paste(x, collapse = " | "))
}

.dv_missing <- function(x) {
  is.na(x) | tolower(trimws(x)) %in% c("", ".", "-", "na", "n/a", "nan", "null", "none")
}

.dv_chrom <- function(x) {
  y <- toupper(sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE))
  y <- sub("^0+(?=[1-9])", "", y, perl = TRUE)
  y[y %in% c("M", "MT", "25", "26")] <- "MT"
  y[y == "23"] <- "X"
  y[y == "24"] <- "Y"
  y[!y %in% c(as.character(1:22), "X", "Y", "MT")] <- NA_character_
  y
}

.dv_contig_to_chrom <- function(contig) {
  out <- rep(NA_character_, length(contig))
  refseq <- grepl("^NC_0000[0-9]{2}\\.", contig)
  num <- as.integer(substr(contig[refseq], 8, 9))
  out[refseq] <- ifelse(num == 23, "X", ifelse(num == 24, "Y", as.character(num)))
  out[grepl("^NC_012920\\.", contig)] <- "MT"
  other <- is.na(out)
  out[other] <- .dv_chrom(contig[other])
  out
}

.dv_assembly <- function(x) {
  b <- tolower(trimws(as.character(x)))
  if (length(b) == 1L && b %in% c("grch38", "hg38")) return("GRCh38")
  if (length(b) == 1L && b %in% c("grch37", "hg19")) return("GRCh37")
  stop("Informe assembly = 'GRCh38' ou 'GRCh37', conforme a montagem da tabela.", call. = FALSE)
}

.dv_term <- function(x) {
  x <- as.character(x)
  if (!length(x)) return(character(0))
  unk <- Encoding(x) == "unknown" & validUTF8(x) & grepl("[^ -~]", x, useBytes = TRUE)
  Encoding(x)[unk] <- "UTF-8"
  x <- chartr("\u00e1\u00e0\u00e2\u00e3\u00e4\u00e9\u00e8\u00ea\u00eb\u00ed\u00ec\u00ee\u00ef\u00f3\u00f2\u00f4\u00f5\u00f6\u00fa\u00f9\u00fb\u00fc\u00e7\u00c1\u00c0\u00c2\u00c3\u00c4\u00c9\u00c8\u00ca\u00cb\u00cd\u00cc\u00ce\u00cf\u00d3\u00d2\u00d4\u00d5\u00d6\u00da\u00d9\u00db\u00dc\u00c7",
              "aaaaaeeeeiiiiooooouuuucAAAAAEEEEIIIIOOOOOUUUUC", trimws(x))
  x <- gsub("[^a-z0-9]+", "_", tolower(x))
  gsub("^_+|_+$", "", x)
}

# "GENE1;GENE2", "a&b", "GENE(dist=10)" -> terms
.dv_split <- function(x) {
  if (length(x) != 1L || .dv_missing(x)) return(character(0))
  parts <- trimws(strsplit(gsub("\\([^)]*\\)", "", x), "[,;|&/]+")[[1L]])
  parts[nzchar(parts) & !.dv_missing(parts)]
}

.dv_find_col <- function(tabela, role, forced) {
  cols <- names(tabela)
  if (!is.null(forced[[role]])) {
    hit <- cols[tolower(cols) == tolower(forced[[role]])]
    if (length(hit) != 1L) stop("Coluna '", forced[[role]], "' (", role, ") nao encontrada.",
                                call. = FALSE)
    return(hit)
  }
  for (a in .dv_ALIASES[[role]]) {
    hit <- cols[tolower(cols) == a]
    if (length(hit)) return(hit[[1L]])
  }
  NULL
}

.dv_vocabulary <- function(extra = NULL) {
  v <- utils::read.delim(text = paste(c("termo\tflags", .dv_VOCAB_TEXT), collapse = "\n"),
                         colClasses = "character", quote = "", stringsAsFactors = FALSE)
  if (!is.null(extra)) {
    e <- utils::read.delim(extra, colClasses = "character", comment.char = "#", quote = "",
                           stringsAsFactors = FALSE)
    names(e)[1:2] <- c("termo", "flags")
    v <- rbind(v[!.dv_term(v$termo) %in% .dv_term(e$termo), ], e[, c("termo", "flags")])
  }
  flags <- lapply(strsplit(v$flags, ",", fixed = TRUE), trimws)
  bad <- setdiff(unlist(flags), c(names(.dv_FLAGS), "INTERGENIC", .dv_UNCHECKABLE))
  if (length(bad)) stop("Flags desconhecidas no vocabulario: ", paste(bad, collapse = ", "),
                        call. = FALSE)
  terms <- .dv_term(v$termo)
  if (anyDuplicated(terms)) stop("Termos repetidos no vocabulario: ",
                                 paste(unique(terms[duplicated(terms)]), collapse = ", "), call. = FALSE)
  stats::setNames(flags, terms)
}

# ---- alleles -------------------------------------------------------------------

.dv_trim_alleles <- function(ref, alt) {
  clean <- function(a) { a <- toupper(trimws(a)); if (a %in% c("-", ".", "")) "" else a }
  r <- clean(ref); a <- clean(alt)
  while (nchar(r) && nchar(a) && substring(r, nchar(r)) == substring(a, nchar(a))) {
    r <- substr(r, 1, nchar(r) - 1); a <- substr(a, 1, nchar(a) - 1)
  }
  while (nchar(r) && nchar(a) && substr(r, 1, 1) == substr(a, 1, 1)) {
    r <- substring(r, 2); a <- substring(a, 2)
  }
  paste0(r, ">", a)
}

.dv_revcomp <- function(s) {
  paste(rev(strsplit(chartr("ACGTacgt", "TGCAtgca", s), "")[[1L]]), collapse = "")
}

# list(ok = TRUE/FALSE/NA, type, note)
.dv_compare_alleles <- function(ref_in, alt_in, db_ref, db_alt) {
  db_pairs <- vapply(strsplit(db_alt, ",", fixed = TRUE)[[1L]],
                     function(a) .dv_trim_alleles(db_ref, a), character(1L))
  db_txt <- paste0(db_ref, ">", db_alt)
  ip <- .dv_trim_alleles(ref_in, alt_in)
  if (ip %in% db_pairs) return(list(ok = TRUE, type = NA_character_, note = NA_character_))
  if (.dv_trim_alleles(.dv_revcomp(toupper(ref_in)), .dv_revcomp(toupper(alt_in))) %in% db_pairs) {
    return(list(ok = FALSE, type = "ALLELES_REVERSE_STRAND",
                note = sprintf("Alelos %s>%s correspondem a fita reversa do dbSNP (%s).",
                                           ref_in, alt_in, db_txt)))
  }
  if (sub(">.*", "", ip) %in% sub(">.*", "", db_pairs)) {
    return(list(ok = NA, type = "ALT_NOT_IN_DBSNP",
                note = sprintf("REF confere, mas ALT %s nao esta entre os alelos do dbSNP (%s).",
                                        alt_in, db_txt)))
  }
  list(ok = FALSE, type = "ALLELES_MISMATCH",
       note = sprintf("Alelos %s>%s nao conferem com o dbSNP (%s).", ref_in, alt_in, db_txt))
}

# ---- dbSNP access ----------------------------------------------------------------

.dv_open <- function(path, assembly, backend = NULL) {
  if (!file.exists(path)) stop("VCF do dbSNP nao encontrado: ", path, call. = FALSE)
  if (!file.exists(paste0(path, ".tbi")) && !file.exists(paste0(path, ".csi"))) {
    stop("Indice nao encontrado (", path, ".tbi). Gere com: tabix -p vcf ", path, call. = FALSE)
  }
  has_rsam <- requireNamespace("Rsamtools", quietly = TRUE)
  has_tabix <- nzchar(Sys.which("tabix"))
  if (is.null(backend)) backend <- if (has_rsam) "Rsamtools" else if (has_tabix) "tabix" else
    stop("Instale Rsamtools (BiocManager::install('Rsamtools')) ou tabix (htslib).", call. = FALSE)
  if (backend == "Rsamtools") {
    tf <- Rsamtools::TabixFile(path)
    header <- Rsamtools::headerTabix(tf)$header
    contigs <- Rsamtools::seqnamesTabix(tf)
  } else {
    header <- system2("tabix", c("-H", shQuote(path)), stdout = TRUE)
    contigs <- system2("tabix", c("-l", shQuote(path)), stdout = TRUE)
  }
  ref_line <- grep("^##reference=", header, value = TRUE)
  vcf_assembly <- if (length(ref_line) && grepl("GRCh38|hg38", ref_line[1], ignore.case = TRUE)) "GRCh38" else
    if (length(ref_line) && grepl("GRCh37|hg19", ref_line[1], ignore.case = TRUE)) "GRCh37" else
    if ("NC_000001.11" %in% contigs) "GRCh38" else if ("NC_000001.10" %in% contigs) "GRCh37" else NA
  if (!is.na(vcf_assembly) && vcf_assembly != assembly) {
    stop(sprintf("O VCF do dbSNP e %s, mas a tabela foi declarada como %s.", vcf_assembly, assembly),
         call. = FALSE)
  }
  chrom <- .dv_contig_to_chrom(contigs)
  keep <- !is.na(chrom) & !duplicated(chrom)
  info_ids <- sub("^##INFO=<ID=([^,]+),.*", "\\1", grep("^##INFO=<ID=", header, value = TRUE))
  build <- sub("^##dbSNP_BUILD_ID=", "", grep("^##dbSNP_BUILD_ID=", header, value = TRUE))
  list(path = path, backend = backend, assembly = vcf_assembly,
       dbsnp_build = if (length(build)) build[[1L]] else NA_character_,
       contig = stats::setNames(contigs[keep], chrom[keep]),
       missing_info = setdiff(.dv_EXPECTED_INFO, info_ids))
}

.dv_parse_vcf <- function(lines) {
  lines <- unique(lines[nzchar(lines) & !startsWith(lines, "#")])
  if (!length(lines)) {
    return(data.frame(chrom = character(0), pos = integer(0), id = character(0),
                      ref = character(0), alt = character(0), vc = character(0),
                      genes = character(0), flags = character(0), stringsAsFactors = FALSE))
  }
  f <- strsplit(lines, "\t", fixed = TRUE)
  get <- function(i) vapply(f, `[`, character(1L), i)
  info <- get(8L)
  value <- function(key) {
    out <- rep("", length(info))
    has <- grepl(paste0("(^|;)", key, "="), info)
    out[has] <- sub(paste0("^;?", key, "="), "",
                    regmatches(info, regexpr(paste0("(^|;)", key, "=[^;]*"), info)))
    out
  }
  symbols <- function(x) vapply(strsplit(x, "|", fixed = TRUE), function(g)
    paste(sub(":.*$", "", g[nzchar(g)]), collapse = ","), character(1L))
  genes <- symbols(value("GENEINFO")); pseudo <- symbols(value("PSEUDOGENEINFO"))
  genes <- ifelse(nzchar(pseudo), ifelse(nzchar(genes), paste(genes, pseudo, sep = ","), pseudo), genes)
  flags <- vapply(strsplit(info, ";", fixed = TRUE), function(t)
    paste(intersect(names(.dv_FLAGS), t), collapse = ","), character(1L))
  data.frame(chrom = .dv_contig_to_chrom(get(1L)), pos = as.integer(get(2L)), id = tolower(get(3L)),
             ref = get(4L), alt = get(5L), vc = value("VC"), genes = genes, flags = flags,
             stringsAsFactors = FALSE)
}

.dv_fetch <- function(ref, chrom, start, end) {
  reg <- unique(data.frame(contig = unname(ref$contig[chrom]), start = as.integer(start),
                           end = as.integer(end), stringsAsFactors = FALSE))
  if (!nrow(reg)) return(character(0))
  if (ref$backend == "Rsamtools") {
    gr <- GenomicRanges::reduce(GenomicRanges::GRanges(reg$contig, IRanges::IRanges(reg$start, reg$end)))
    lines <- unlist(Rsamtools::scanTabix(Rsamtools::TabixFile(ref$path), param = gr), use.names = FALSE)
  } else {
    tmp <- tempfile(fileext = ".txt"); on.exit(unlink(tmp), add = TRUE)
    utils::write.table(reg[order(reg$contig, reg$start), ], tmp, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)
    lines <- system2("tabix", c("-R", shQuote(tmp), shQuote(ref$path)), stdout = TRUE)
  }
  unique(lines)
}

# Reads the WHOLE VCF (tens of minutes on full dbSNP); only for rsIDs not found at POS.
.dv_scan <- function(ref, rsids) {
  rsids <- unique(rsids[!is.na(rsids)])
  if (!length(rsids)) return(.dv_parse_vcf(character(0)))
  if (!nzchar(Sys.which("awk"))) stop("--scan-rsids precisa de awk (Linux/macOS/WSL).", call. = FALSE)
  ids <- tempfile(fileext = ".txt"); on.exit(unlink(ids), add = TRUE)
  writeLines(rsids, ids)
  decomp <- if (nzchar(Sys.which("bgzip"))) "bgzip -dc -@ 4" else "gzip -dc"
  cmd <- sprintf("%s %s | awk -F'\\t' 'NR==FNR{ids[$1];next} !/^#/ && (tolower($3) in ids)' %s -",
                 decomp, shQuote(ref$path), shQuote(ids))
  .dv_parse_vcf(system(cmd, intern = TRUE))
}

# ---- functional terms ------------------------------------------------------------

.dv_describe <- function(flags, genes) {
  if (!length(flags)) return(if (!length(genes)) "intergenic_variant" else "(sem anotacao funcional)")
  paste(unname(.dv_FLAGS[flags]), collapse = ", ")
}

# list(check = TRUE/FALSE/NA, notes = list(type, level, msg))  (character vectors)
.dv_check_terms <- function(terms, label, flags, genes, vocab) {
  notes <- list(type = character(0), level = character(0), msg = character(0))
  add <- function(type, level, msg) {
    notes$type  <<- c(notes$type, rep_len(paste0(label, "_", type), length(msg)))
    notes$level <<- c(notes$level, rep_len(level, length(msg)))
    notes$msg   <<- c(notes$msg, msg)
  }
  if (!length(terms)) return(list(check = NA, notes = notes))
  db <- if (!length(flags) && !length(genes)) "INTERGENIC" else flags
  if (any(c("ASS", "DSS") %in% db)) db <- union(db, "INT")  # splice sites lie in the intron
  norm <- .dv_term(terms)
  known <- norm %in% names(vocab)
  checkable <- known & !vapply(norm, function(t) identical(vocab[[t]], .dv_UNCHECKABLE), logical(1L))
  if (any(!known)) add("UNKNOWN_TERM", "REVIEW",
                       sprintf("%s: termo '%s' fora do vocabulario.", label, terms[!known]))
  if (!any(checkable)) return(list(check = NA, notes = notes))
  if (!length(db)) {
    add("NOT_VERIFIABLE", "REVIEW", sprintf("%s: dbSNP sem anotacao funcional; nao conferido.", label))
    return(list(check = NA, notes = notes))
  }
  ok <- vapply(norm[checkable], function(t) any(vocab[[t]] %in% db), logical(1L))
  if (any(!ok)) add("MISMATCH", "FAIL", sprintf("%s: '%s' nao confere com o dbSNP (%s).", label,
                                                terms[checkable][!ok], .dv_describe(flags, genes)))
  list(check = all(ok), notes = notes)
}

.dv_effect_vs_impact <- function(effect_terms, impact_terms, vocab) {
  if (!length(effect_terms) || !length(impact_terms)) return(NA_character_)
  cls <- unique(stats::na.omit(.dv_CLASS_TERMS[.dv_term(impact_terms)]))
  if (length(cls) != 1L) return(NA_character_)
  flags <- setdiff(unlist(lapply(.dv_term(effect_terms), function(t) vocab[[t]])), .dv_UNCHECKABLE)
  if (!length(flags)) return(NA_character_)
  expected <- unique(unname(.dv_CLASS[flags]))
  if (cls %in% expected) return(NA_character_)
  sprintf("IMPACT %s incoerente com EFFECT '%s' (esperado: %s).", cls,
          paste(effect_terms, collapse = "&"), paste(expected, collapse = " ou "))
}

# ---- input ------------------------------------------------------------------------

.dv_read <- function(path) {
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

# ---- main function ---------------------------------------------------------------

validar_dbsnp <- function(tabela, assembly, dbsnp, colunas = list(), janela = 10L,
                          scan_rsids = FALSE, vocabulario = NULL, backend = NULL) {
  assembly <- .dv_assembly(if (missing(assembly)) NA else assembly)
  if (missing(dbsnp)) stop("Informe dbsnp = caminho do VCF do dbSNP (.gz com .tbi).", call. = FALSE)
  janela <- suppressWarnings(as.integer(janela))
  if (length(janela) != 1L || is.na(janela) || janela < 0L) stop("janela deve ser inteiro >= 0.", call. = FALSE)
  if (is.character(tabela) && length(tabela) == 1L) tabela <- .dv_read(tabela)
  if (!is.data.frame(tabela)) stop("tabela deve ser um data.frame ou caminho de TSV.", call. = FALSE)
  tabela <- as.data.frame(tabela, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(tolower(names(tabela)))) stop("Ha nomes de colunas duplicados na entrada.", call. = FALSE)

  cols <- lapply(stats::setNames(nm = names(.dv_ALIASES)), .dv_find_col, tabela = tabela, forced = colunas)
  if (is.null(cols$chrom) || is.null(cols$pos)) stop("Colunas obrigatorias: CHROM e POS.", call. = FALSE)
  col_value <- function(role) if (is.null(cols[[role]])) rep(NA_character_, nrow(tabela)) else
    as.character(tabela[[cols[[role]]]])

  ref <- .dv_open(dbsnp, assembly, backend)
  vocab <- .dv_vocabulary(vocabulario)
  n <- nrow(tabela)

  added <- c("DBSNP_STATUS", "DBSNP_RSID", "DBSNP_POS", "DBSNP_REF", "DBSNP_ALT", "DBSNP_VC",
             "DBSNP_GENES", "DBSNP_REGION", "DBSNP_CONSEQUENCE", "CHECK_RSID", "CHECK_ALLELES",
             "CHECK_GENE", "CHECK_REGION", "CHECK_EFFECT", "CHECK_IMPACT", "DBSNP_VALIDATION",
             "DBSNP_NOTES")
  collision <- intersect(tolower(names(tabela)), tolower(added))
  if (length(collision)) stop("A entrada ja contem colunas de validacao: ", paste(collision, collapse = ", "),
                              call. = FALSE)
  # plain vectors while filling (a data.frame cell-by-cell is very slow in R)
  out <- lapply(stats::setNames(nm = added), function(k)
    if (startsWith(k, "CHECK_")) rep(NA, n) else rep(NA_character_, n))
  # every note is logged with a type and a level (FAIL, REVIEW or INFO)
  log_chunks <- vector("list", 1024L); n_chunks <- 0L
  add_note <- function(i, type, level, msg) {
    k <- max(length(i), length(msg))
    n_chunks <<- n_chunks + 1L
    if (n_chunks > length(log_chunks)) length(log_chunks) <<- 2L * length(log_chunks)
    log_chunks[[n_chunks]] <<- list(rep_len(i, k), rep_len(type, k), rep_len(level, k), rep_len(msg, k))
  }

  # normalised keys used only for comparison
  chrom <- .dv_chrom(col_value("chrom"))
  pos_txt <- trimws(col_value("pos"))
  pos <- ifelse(grepl("^[0-9]{1,9}$", pos_txt), suppressWarnings(as.integer(pos_txt)), NA_integer_)
  rs <- tolower(trimws(col_value("rsid")))
  rs[!grepl("^rs[1-9][0-9]*$", rs)] <- NA_character_
  ref_in <- col_value("ref"); alt_in <- col_value("alt")
  allele_missing <- function(x) is.na(x) | (.dv_missing(x) & trimws(x) != "-")
  has_alleles <- !is.null(cols$ref) && !is.null(cols$alt)
  alleles_ok <- has_alleles & !allele_missing(ref_in) & !allele_missing(alt_in) &
    !(trimws(ref_in) %in% "-" & trimws(alt_in) %in% "-")
  split_col <- function(role) {
    v <- col_value(role); u <- unique(v)
    parsed <- lapply(u, .dv_split)
    function(i) parsed[[match(v[[i]], u)]]
  }
  terms <- lapply(c(gene = "gene", region = "region", effect = "effect", impact = "impact"), split_col)

  # one indexed query for the whole table
  ok_site <- !is.na(chrom) & !is.na(pos) & chrom %in% names(ref$contig)
  rec <- .dv_parse_vcf(.dv_fetch(ref, chrom[ok_site], pmax(1L, pos[ok_site] - janela), pos[ok_site] + janela))
  rec <- rec[order(rec$chrom, rec$pos), , drop = FALSE]
  by_chrom <- split(seq_len(nrow(rec)), rec$chrom)
  pos_by_chrom <- lapply(by_chrom, function(k) as.numeric(rec$pos[k]))
  near_of <- function(c, p) {
    idx <- by_chrom[[c]]; if (is.null(idx)) return(integer(0))
    q <- pos_by_chrom[[c]]
    a <- findInterval(p - janela - 1, q) + 1L; b <- findInterval(p + janela, q)
    if (b >= a) idx[a:b] else integer(0)
  }

  chosen_of <- rep(NA_integer_, n)
  to_scan <- integer(0)
  gene_bad <- logical(n)

  for (i in seq_len(n)) {
    if (is.na(chrom[i]) || is.na(pos[i])) {
      out$DBSNP_STATUS[i] <- "NOT_EVALUATED"
      add_note(i, "POSITION_UNREADABLE", "INFO", "CHROM/POS nao legiveis; linha nao comparada.")
      next
    }
    if (!chrom[i] %in% names(ref$contig)) {
      out$DBSNP_STATUS[i] <- "CONTIG_NOT_IN_DBSNP"
      add_note(i, "CONTIG_NOT_IN_DBSNP", "INFO", sprintf("Cromossomo %s ausente no VCF do dbSNP usado.", chrom[i]))
      next
    }
    near <- near_of(chrom[i], pos[i])
    exact <- near[rec$pos[near] == pos[i]]
    pick <- function(cands) {
      if (length(cands) <= 1L || !alleles_ok[i]) return(cands[1L])
      hit <- vapply(cands, function(k) isTRUE(.dv_compare_alleles(ref_in[i], alt_in[i], rec$ref[k], rec$alt[k])$ok),
                    logical(1L))
      if (any(hit)) cands[hit][1L] else cands[1L]
    }
    chosen <- NA_integer_
    if (!is.na(rs[i])) {
      same <- near[rec$id[near] == rs[i]]
      if (length(same)) {
        chosen <- same[which.min(abs(rec$pos[same] - pos[i]))]
        out$CHECK_RSID[i] <- TRUE
        if (rec$pos[chosen] == pos[i]) {
          out$DBSNP_STATUS[i] <- "CONFIRMED"
        } else {
          out$DBSNP_STATUS[i] <- "NEIGHBOR_POSITION"
          add_note(i, "NEIGHBOR_POSITION", "REVIEW",
                   sprintf("%s esta em %s:%d no dbSNP (%+d bp); convencao de indel diferente?",
                           rs[i], chrom[i], rec$pos[chosen], rec$pos[chosen] - pos[i]))
        }
      } else if (length(exact)) {
        chosen <- pick(exact)
        out$DBSNP_STATUS[i] <- "RSID_MISMATCH"
        out$CHECK_RSID[i] <- FALSE
        add_note(i, "RSID_MISMATCH", "FAIL", sprintf("No dbSNP %s:%d corresponde a %s, nao a %s.",
                 chrom[i], pos[i], paste(rec$id[exact], collapse = ","), rs[i]))
      } else {
        out$DBSNP_STATUS[i] <- "NOT_AT_POSITION"
        out$CHECK_RSID[i] <- FALSE
        to_scan <- c(to_scan, i)
      }
    } else if (length(exact)) {
      chosen <- pick(exact)
      out$DBSNP_STATUS[i] <- "RSID_AVAILABLE"
      if (!is.null(cols$rsid)) {
        add_note(i, "RSID_AVAILABLE", "REVIEW", sprintf("Linha sem rsID; o dbSNP tem %s nesta posicao.",
                                                        paste(rec$id[exact], collapse = ",")))
      }
    } else {
      out$DBSNP_STATUS[i] <- "NOVEL"
    }
    if (is.na(chosen)) next
    chosen_of[i] <- chosen

    genes <- if (nzchar(rec$genes[chosen])) strsplit(rec$genes[chosen], ",")[[1L]] else character(0)
    flags <- if (nzchar(rec$flags[chosen])) strsplit(rec$flags[chosen], ",")[[1L]] else character(0)
    out$DBSNP_RSID[i] <- rec$id[chosen]; out$DBSNP_POS[i] <- as.character(rec$pos[chosen])
    out$DBSNP_REF[i] <- rec$ref[chosen]; out$DBSNP_ALT[i] <- rec$alt[chosen]
    out$DBSNP_VC[i] <- rec$vc[chosen]; out$DBSNP_GENES[i] <- rec$genes[chosen]
    out$DBSNP_REGION[i] <- if (length(flags)) paste(unique(unname(.dv_REGION[flags])), collapse = ",") else
      if (!length(genes)) "intergenic" else NA_character_
    out$DBSNP_CONSEQUENCE[i] <- if (length(flags)) paste(unname(.dv_FLAGS[flags]), collapse = ",") else
      if (!length(genes)) "intergenic_variant" else NA_character_

    if (alleles_ok[i]) {
      cmp <- .dv_compare_alleles(ref_in[i], alt_in[i], rec$ref[chosen], rec$alt[chosen])
      out$CHECK_ALLELES[i] <- cmp$ok
      if (!is.na(cmp$note)) add_note(i, cmp$type, if (isFALSE(cmp$ok)) "FAIL" else "REVIEW", cmp$note)
    }

    if (!is.null(cols$gene)) {
      g <- terms$gene(i)
      if (!length(g)) {
        if (length(genes)) {
          add_note(i, "GENE_EMPTY", "REVIEW", sprintf("GENE vazio; dbSNP associa %s.", paste(genes, collapse = ",")))
        }
      } else {
        bad <- g[!toupper(g) %in% toupper(genes)]
        out$CHECK_GENE[i] <- !length(bad)
        if (length(bad)) {
          gene_bad[i] <- TRUE
          add_note(i, "GENE_MISMATCH", "FAIL", if (length(genes))
            sprintf("GENE %s nao confere com o dbSNP (%s); simbolo antigo/alias?",
                    paste(bad, collapse = ","), paste(genes, collapse = ",")) else
            sprintf("dbSNP nao associa gene a esta variante, mas a linha indica %s.", paste(bad, collapse = ",")))
        }
      }
    }

    for (role in c("region", "effect", "impact")) {
      if (is.null(cols[[role]])) next
      r <- .dv_check_terms(terms[[role]](i), toupper(role), flags, genes, vocab)
      out[[paste0("CHECK_", toupper(role))]][i] <- r$check
      if (length(r$notes$msg)) add_note(i, r$notes$type, r$notes$level, r$notes$msg)
    }
  }

  # rows of the same dbSNP record (one row per gene): dbSNP gene without a row
  if (!is.null(cols$gene)) {
    linked <- which(!is.na(chosen_of) & nzchar(rec$genes[pmax(chosen_of, 1L)]))
    for (g in split(linked, chosen_of[linked])) {
      if (any(gene_bad[g])) next
      in_table <- unique(toupper(unlist(lapply(g, terms$gene))))
      if (!length(in_table)) next
      db_genes <- strsplit(rec$genes[chosen_of[g[1L]]], ",")[[1L]]
      absent <- db_genes[!toupper(db_genes) %in% in_table]
      if (length(absent)) {
        add_note(g, "GENE_ROW_MISSING", "REVIEW", sprintf("dbSNP tambem associa %s a %s; nao ha linha para esse gene.",
                                                          rec$id[chosen_of[g[1L]]], paste(absent, collapse = ",")))
      }
    }
  }

  # IMPACT class vs EFFECT in the same row
  if (!is.null(cols$effect) && !is.null(cols$impact)) {
    for (i in seq_len(n)) {
      msg <- .dv_effect_vs_impact(terms$effect(i), terms$impact(i), vocab)
      if (!is.na(msg)) add_note(i, "EFFECT_IMPACT_INCONSISTENT", "REVIEW", msg)
    }
  }

  # rsIDs not found at POS
  if (length(to_scan)) {
    found <- if (scan_rsids) .dv_scan(ref, rs[to_scan]) else NULL
    for (i in to_scan) {
      if (is.null(found)) {
        add_note(i, "RSID_NOT_AT_POSITION", "FAIL",
                 sprintf("%s nao encontrado em %s:%d (+-%d bp); use --scan-rsids para procurar no arquivo todo.",
                         rs[i], chrom[i], pos[i], janela))
      } else if (rs[i] %in% found$id) {
        k <- which(found$id == rs[i])[1L]
        out$DBSNP_STATUS[i] <- "RSID_ELSEWHERE"
        add_note(i, "RSID_ELSEWHERE", "FAIL", sprintf("%s esta em %s:%d (%s>%s) no dbSNP.", rs[i],
                 found$chrom[k], found$pos[k], found$ref[k], found$alt[k]))
      } else {
        out$DBSNP_STATUS[i] <- "RSID_NOT_IN_DBSNP"
        add_note(i, "RSID_NOT_IN_DBSNP", "FAIL", sprintf("%s nao existe neste dbSNP (build %s); mesclado ou retirado?",
                                                         rs[i], ref$dbsnp_build))
      }
    }
  }

  out <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
  checks <- as.matrix(out[grep("^CHECK_", names(out))])
  any_false <- apply(checks, 1L, function(r) any(r %in% FALSE))
  any_true <- apply(checks, 1L, function(r) any(r %in% TRUE))
  chunks <- log_chunks[seq_len(n_chunks)]
  log_row   <- as.integer(unlist(lapply(chunks, `[[`, 1L)))
  log_type  <- as.character(unlist(lapply(chunks, `[[`, 2L)))
  log_level <- as.character(unlist(lapply(chunks, `[[`, 3L)))
  log_msg   <- as.character(unlist(lapply(chunks, `[[`, 4L)))
  review <- seq_len(n) %in% log_row[log_level == "REVIEW"]
  out$DBSNP_VALIDATION <- ifelse(any_false, "FAIL", ifelse(review, "REVIEW",
                           ifelse(any_true, "PASS", "NOT_EVALUATED")))
  by_row <- split(log_msg, factor(log_row, levels = seq_len(n)))
  out$DBSNP_NOTES <- vapply(by_row, .dv_note, character(1L), USE.NAMES = FALSE)

  # one row per problem (FAIL/REVIEW), for a separate report
  keep <- log_level %in% c("FAIL", "REVIEW")
  ord <- order(log_row[keep], match(log_level[keep], c("FAIL", "REVIEW")))
  r <- log_row[keep][ord]
  problems <- data.frame(ROW = r, CHROM = col_value("chrom")[r], POS = col_value("pos")[r],
                         ID = col_value("rsid")[r], REF = col_value("ref")[r],
                         ALT = col_value("alt")[r], GENE = col_value("gene")[r],
                         DBSNP_VALIDATION = out$DBSNP_VALIDATION[r], LEVEL = log_level[keep][ord],
                         TYPE = log_type[keep][ord], MESSAGE = gsub("[\t\r\n]+", " ", log_msg[keep][ord]),
                         DBSNP_RSID = out$DBSNP_RSID[r], DBSNP_POS = out$DBSNP_POS[r],
                         stringsAsFactors = FALSE, check.names = FALSE)

  result <- cbind(tabela, out)
  with_rs <- !is.na(out$CHECK_RSID)
  if (sum(with_rs) >= 10L && mean(!out$CHECK_RSID[with_rs]) > 0.5) {
    warning(sprintf("%.0f%% dos rsIDs nao conferem com a posicao no dbSNP %s. A tabela esta em outro build?",
                    100 * mean(!out$CHECK_RSID[with_rs]), ref$assembly), call. = FALSE)
  }
  attr(result, "dbsnp_file") <- ref$path
  attr(result, "dbsnp_build") <- ref$dbsnp_build
  attr(result, "assembly") <- assembly
  attr(result, "columns_used") <- unlist(cols)
  attr(result, "dbsnp_missing_info") <- ref$missing_info
  attr(result, "problems") <- problems
  result
}

salvar_dbsnp_tsv <- function(tabela, arquivo) {
  if (file.exists(arquivo)) stop("A saida ja existe: ", arquivo, call. = FALSE)
  if (!dir.exists(dirname(arquivo))) stop("A pasta de saida nao existe.", call. = FALSE)
  if (any(grepl("[\t\r\n]", names(tabela))) ||
      any(vapply(tabela, function(x) any(grepl("[\t\r\n]", as.character(x))), logical(1L)))) {
    stop("Ha campos com tabulacao/quebra de linha; revise antes de exportar TSV.", call. = FALSE)
  }
  utils::write.table(tabela, file = arquivo, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
  invisible(arquivo)
}

# One row per problem (FAIL/REVIEW) of a validar_dbsnp() result.
# ROW is the data row of the input table (1 = first row after the header).
problemas_dbsnp <- function(resultado) {
  p <- attr(resultado, "problems")
  if (is.null(p)) stop("Use o resultado de validar_dbsnp().", call. = FALSE)
  p
}

# ---- command line ---------------------------------------------------------------

.dv_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste("Uso: Rscript validar_dbsnp.R entrada.tsv saida.tsv GRCh38|GRCh37 --dbsnp ARQ.gz",
                 "       [--problemas ARQ.tsv] [--scan-rsids] [--janela N] [--vocabulario ARQ.tsv]",
                 "       [--col papel=COLUNA ...]",
                 sep = "\n")
  if (any(tolower(args) %in% c("-h", "--help"))) { cat(usage, "\n"); return(invisible(NULL)) }
  opt <- list(dbsnp = NULL, janela = 10L, vocabulario = NULL, problemas = NULL, scan = FALSE,
              colunas = list())
  pos_args <- character(0)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]; key <- tolower(a)
    if (key == "--scan-rsids") {
      opt$scan <- TRUE
    } else if (key %in% c("--dbsnp", "--janela", "--vocabulario", "--problemas", "--col")) {
      if (i == length(args)) stop("Falta o valor de ", a, call. = FALSE)
      v <- args[[i + 1L]]; i <- i + 1L
      if (key == "--col") {
        kv <- strsplit(v, "=", fixed = TRUE)[[1L]]
        role <- tolower(kv[[1L]])
        if (length(kv) != 2L || !role %in% names(.dv_ALIASES)) {
          stop("--col espera papel=COLUNA, papel em: ", paste(names(.dv_ALIASES), collapse = ", "), call. = FALSE)
        }
        opt$colunas[[role]] <- kv[[2L]]
      } else opt[[sub("^--", "", key)]] <- v
    } else if (startsWith(a, "--")) {
      stop("Opcao desconhecida: ", a, "\n", usage, call. = FALSE)
    } else pos_args <- c(pos_args, a)
    i <- i + 1L
  }
  if (length(pos_args) != 3L || is.null(opt$dbsnp)) stop(usage, call. = FALSE)
  if (is.null(opt$problemas)) opt$problemas <- paste0(sub("\\.tsv$", "", pos_args[[2L]]), "_problemas.tsv")
  for (f in c(pos_args[[2L]], opt$problemas)) {
    if (file.exists(f)) stop("A saida ja existe: ", f, call. = FALSE)
    if (!dir.exists(dirname(f))) stop("A pasta de saida nao existe: ", dirname(f), call. = FALSE)
  }
  if (identical(normalizePath(pos_args[[2L]], mustWork = FALSE), normalizePath(opt$problemas, mustWork = FALSE))) {
    stop("--problemas precisa ser diferente da saida principal.", call. = FALSE)
  }
  result <- validar_dbsnp(pos_args[[1L]], assembly = pos_args[[3L]], dbsnp = opt$dbsnp,
                          colunas = opt$colunas, janela = opt$janela, scan_rsids = opt$scan,
                          vocabulario = opt$vocabulario)
  salvar_dbsnp_tsv(result, pos_args[[2L]])
  problems <- problemas_dbsnp(result)
  salvar_dbsnp_tsv(problems, opt$problemas)
  used <- attr(result, "columns_used")
  message("dbSNP build ", attr(result, "dbsnp_build"), " | colunas: ",
          paste(names(used), used, sep = "=", collapse = ", "))
  if (length(attr(result, "dbsnp_missing_info"))) {
    message("ATENCAO: VCF sem os campos INFO ", paste(attr(result, "dbsnp_missing_info"), collapse = ", "))
  }
  print(table(result$DBSNP_STATUS, useNA = "ifany"))
  print(table(result$DBSNP_VALIDATION, useNA = "ifany"))
  if (nrow(problems)) {
    message("Problemas por tipo:")
    tt <- table(paste(problems$LEVEL, problems$TYPE))
    print(tt[order(-as.integer(tt))])
  }
  message("Tabela salva: ", pos_args[[2L]])
  message("Problemas (", nrow(problems), " em ", length(unique(problems$ROW)), " linhas): ", opt$problemas)
  invisible(result)
}

if (sys.nframe() == 0L) {
  tryCatch(.dv_main(), error = function(e) {
    message("ERRO: ", conditionMessage(e))
    quit(save = "no", status = 1L)
  })
}

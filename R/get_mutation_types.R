#' Gets a list of the standard 96 mutation types used in SNV mutation signature analysis.
#'
#' @return A 96-element character vector containing SNV mutation types
#'
#' @import dplyr
#' @import tidyr
#' @export

all_snv_mutation_types <- function() {
    bases = c('A', 'C', 'G', 'T')
    changes = c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
    all_mutation_types <- crossing(
        five_prime = bases,
        substitution = changes,
        three_prime = bases
    ) %>% 
    mutate(
        mutation_type = paste0(five_prime, '[', substitution, ']', three_prime)
    ) %>% 
    .$mutation_type
}

#' Gets the trinucleotide context for a set of genomic coordinates.
#'
#' @param chr Character vector of chromosome names. Must match provided genome.
#' @param pos Integer vector of genomic positions.
#' @param genome BSgenome object - default is BSgenome.Hsapiens.UCSC.hg19
#' @return Character vector of trinucleotide contexts
#'
#' @import BSgenome
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import dplyr

get_trinucleotide <- function(chr, pos, ref, genome = BSgenome.Hsapiens.UCSC.hg19) {
    chr_levels <- chr %>% unique
    if (all(chr_levels %in% seqlevels(genome))) {
        message('Chromosome levels match. Proceeding...')
        chromosome = chr
    } else if (all(paste0('chr', chr_levels) %in% seqlevels(genome))) {
        message('Chromosome levels match after prepending "chr". Proceeding...')
        warning('"chr" was automatically prepended to all chromosome names to match genome seqlevels')
        chromosome = paste0('chr', chr)
    } else if (all(gsub('chr', '', chr_levels) %in% seqlevels(genome))) {
        message('Chromosome levels match after removing "chr". Proceeding...')
        warning('"chr" was automatically removed from all chromosome names to match genome seqlevels')
        chromosome = gsub('chr', '', chr)
    } else {
        stop('Some chromosome names do not match genome seqlevels. Cannot proceed.')
    }

    trinucleotide = getSeq(genome, chromosome, pos - 1, pos + 1) %>% as.character
    genome_ref = substr(trinucleotide, 2, 2)
    
    if (any(genome_ref != ref)) {
        if (all(genome_ref[genome_ref != ref] == 'N')) {
            warning(sprintf('%s positions matched N in genome ref.', sum(genome_ref != ref)))
        } else {
            stop('Some positions did not match reference genome. Stopping.')
        }
    }

    return(trinucleotide)
}

#' Returns the complementary sequence for a set of strings
#'
#' @param base_strings  Character vector of sequences to return the complement of.
#'                      Any characters aside from A, C, G, T will be returned as is.
#'
#' @param reverse       If TRUE, then the sequence is reversed before it is returned.
#'
#' @return Character vector of complement or reverse complement sequences
#'
#' @import dplyr

base_complement <- function(base_strings, reverse = FALSE) {
    conversion = c(A = 'T', T = 'A', C = 'G', G = 'C')
    strsplit(base_strings, '') %>%
        sapply(function(bases) {
            base_vec = if_else(
                bases %in% names(conversion),
                conversion[bases],
                bases
            )
            if (reverse) {
                base_vec <- rev(base_vec)
            }
            base_vec %>% paste(collapse='') 
        })
}

#' Collapses mutation types of G>N and A>N mutations into C>N and T>N mutations
#' so that there are 6 possible base changes.
#'
#' @param mutation_types Character vector of mutation types in format C[T>A]G
#'                       as an example of a T>A mutation in CTG trinucleotide context.
#' 
#' @return New mutation types with any G>N/A>N mutations replaced and their
#'         trinucleotide contexts reverse complemented
#'

collapse_mutation_types <- function(mutation_types) {
    trinucleotide = gsub('(.)\\[(.)>(.)\\](.)', '\\1\\2\\4', mutation_types)
    substitution = gsub('(.)\\[(.)>(.)\\](.)', '\\2>\\3', mutation_types)
    ref_base = gsub('(.)\\[(.)>(.)\\](.)', '\\2', mutation_types)

    to_reverse = ref_base %in% c('G', 'A')
    trinucleotide[to_reverse] = base_complement(trinucleotide[to_reverse], reverse = TRUE)
    substitution[to_reverse] = base_complement(substitution[to_reverse], reverse = FALSE)

    stitch_mutation_types(substitution, trinucleotide)
}

#' Given a set of substitutions and trinucleotides, merge into unique mutation
#' types. For example, T>A mutation in CTG context will become C[T>A]G.
#'
#' @param substitution Character vector of substitutions in the form T>A
#' @param trinucleotide Character vector of trinucleotide contexts in the form CTG
#' @return Character vector of mutation types in the form C[T>A]G

stitch_mutation_types <- function(substitution, trinucleotide) {
    sprintf(
        '%s[%s]%s',
        substr(trinucleotide, 1, 1),
        substitution,
        substr(trinucleotide, 3, 3)
    )
}

#' Given chromosome, position, reference, and altered allele values, return a vector of
#' SNV mutation types conforming to the standard 96-element mutation types.
#'
#' @param chr       Character vector of chromosomes (must match genome seqlevels)
#' @param pos       Integer vector of mutation positions
#' @param ref       Character vector of reference alleles (A, C, G, or T)
#' @param alt       Character vector of mutated alleles (A, C, G, or T)
#' @param genome    BSgenome object from which to derive trinucleotide context sequences
#' @return          Character vector of mutation types
#'
#' @import tibble
#' @import dplyr
#' @export

get_snv_mutation_type <- function(chr, pos, ref, alt, genome = BSgenome.Hsapiens.UCSC.hg19) {
    tibble(chr, pos, ref, alt) %>%
        mutate(
            trinucleotide_ = get_trinucleotide(chr, pos, ref, genome = genome),
            substitution_ = paste(ref, alt, sep='>'),
            mutation_type_ = stitch_mutation_types(substitution_, trinucleotide_)
        ) %>%
        mutate(mutation_type = collapse_mutation_types(mutation_type_)) %>%
        .$mutation_type
}

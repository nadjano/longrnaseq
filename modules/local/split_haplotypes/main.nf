process SPLIT_HAPLOTYPES {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"
    
    input:
    tuple val(meta), path(input_gff)
    
    output:
    tuple val(meta), path("${meta.id}.*.gff"), emit: split_gffs
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract unique haplotype suffixes from the GFF file.
    haplotypes=\$(grep -o '[Cc]hr[0-9]\\+_\\?[A-Z0-9]' "$input_gff" | \
    sed 's/.*\\([A-Z0-9]\\)\$/\\1/' | \
    sort -u)

    # Split into haplotype-specific files
    for hap in \$haplotypes; do
        echo \$hap
        if [[ "\$hap" =~ ^[0-9]+\$ ]]; then
            # Numeric haplotype (original pattern)
            grep "chr_\\?[0-9]\\+_\${hap}\\b" "$input_gff" > "${meta.id}.hap\${hap}.gff"
        else
            # Alphanumeric haplotype (new pattern)
            grep "[Cc]hr[0-9]\\+_\\?\${hap}\\b" "$input_gff" > "hap\${hap}.gff"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | head -n1 | sed 's/.*[[:space:]]\\([0-9][0-9.]*\\).*/\\1/')
    END_VERSIONS
    """
}

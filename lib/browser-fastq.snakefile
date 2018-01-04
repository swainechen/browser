shell.executable("/bin/bash")
# "unofficial bash strict mode" http://www.redsymbol.net/articles/unofficial-bash-strict-mode/
shell.prefix("set -euo pipefail;")

assert 'ACC' in config, ("required --config ACC=xyz missing")

onerror:
    print("FASTQ PREP ERROR")
    
onsuccess:
    print("FASTQ PREP SUCCESS")

# first rule defines final target
rule cleanup:
    shell:
        "if {config[GETFILES]} -checkonly {config[ACC]}; then rm -f `{config[GETFILES]} {config[ACC]} -delimiter \" \"`; fi;"

rule fastq:
    shell:
        "{config[GETFILES]} {config[ACC]} && {config[IMPORT_FASTQ]} {config[ACC]} -db;"

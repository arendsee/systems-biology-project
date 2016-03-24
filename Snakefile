rule getAdviceFromPonies:
    message:
        "Ponies forever ..."
    shell:
        "for p in `ponysay -l`; do ponysay -of $p; done"

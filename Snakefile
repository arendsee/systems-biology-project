stuid = [
    'SRP009850', # shoot apical meristems stages
    'SRP018034', # leaf senesence stages
    'SRP053394', # hi/lo temp
    'SRP059384', # roots versus shoots
    'SRP059724', # (hi/lo P) X (mock/Ct/Ci) X (6,10,16,24 dpi)
    'SRP063314', # leaf, flower, root
    'SRP064782', # circadian clock ???
    'SRP063421', # (hi/lo Pi) X (root, shoot) X (long, short exposure)
    'SRP069266'  # high altitude adaptation
]

# make a dictionary of runids? relationships?
# do each in its own dedicated script
# then write a merge script to pool like an idiot

rule getAdviceFromPonies:
    message:
        "Ponies forever ..."
    shell:
        "for p in `ponysay -l`; do ponysay -of $p; done"

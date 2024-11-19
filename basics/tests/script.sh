#!/bin/bash

# Liste des dossiers à déplacer
folders=(
    "./survival"
    "./codetools"
)

# Répertoire de destination
destination="/gpfs/workdir/selvestra"

# Déplacer chaque dossier
for folder in "${folders[@]}"; do
    mv "$folder" "$destination"
done
#!/bin/bash

# ===== USER SETTINGS =====
REMOTE_URL="https://github.com/Jojo6297/Metabolic_set_theory.git"
# =========================

# Initialize git repo if not already initialized
if [ ! -d ".git" ]; then
    git init
    git branch -M main
    git remote add origin "$REMOTE_URL"
fi

# Add all files
git add .

# Commit with an empty commit message
git commit --allow-empty-message -m ""

# Push to GitHub
git push -u origin main
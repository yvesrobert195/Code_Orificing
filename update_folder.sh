#!/bin/bash

git add -u
git reset -- Clones/*
git reset -- User_Input.m
git reset -- *~
git reset -- */*~
git reset -- */*/*~

if [ $# -eq 0 ]
then
    git commit -m "Update without comment"
    echo "No arguments supplied"
else
    git commit -m "$1"
fi

git push

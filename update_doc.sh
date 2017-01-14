#!/bin/bash

git fetch origin gh-pages
git checkout gh-pages
if [ $? -eq 0 ]; then
    make generate_doc
    cp -R doc/html/* .
    git add *.html
    git add standalone_HTML.manifest
    git add Aboria/*
    git add aboria/*
    git add Eigen/*
    git add index/*
    git commit -m "Updated documentation by TravisCI"
    git push
    git checkout master
    echo OK
else
    echo FAIL
fi


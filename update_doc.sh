#!/bin/bash

make generate_doc
git fetch origin gh-pages
git checkout gh-pages
if [ $? -eq 0 ]; then
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


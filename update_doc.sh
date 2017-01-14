#!/bin/bash

git fetch origin gh-pages
git checkout gh-pages
make generate_doc
cp -R doc/html/* .
git add *.html
git add Aboria/*
git add aboria/*
git add Eigen/*
git add index/*
git commit -m "Updated documentation by TravisCI"
git push
git checkout master

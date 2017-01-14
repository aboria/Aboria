#!/usr/bash

make generate_doc
git checkout gh-pages
cp -R doc/html/* .
git add *.html
git add Aboria/*
git add aboria/*
git add Eigen/*
git add index/*
git commit -m "Updated documentation by TravisCI"
git push

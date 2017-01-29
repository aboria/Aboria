#!/bin/bash

make aboria-html-website
git fetch origin gh-pages
git checkout gh-pages
if [ $? -eq 0 ]; then
    cp -R doc/html-website/* .
    git add *.html
    git add aboria/\*.html
    git add Aboria/\*.html
    git add css/\*.css
    git add images/\*.png
    git add images/\*.jpg
    git add images/\*.svg
    git add HTML.manifest
    git commit -m "Updated documentation by TravisCI"
    git push
    git checkout master
    echo 'UPDATE_DOC SUCCEEDED :)'
else
    echo 'UPDATE_DOC FAILED :('
fi


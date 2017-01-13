#!/usr/bash

make generate_doc
checkout gh-pages
cp -R doc/html/* .
git add *
git commit -m "Updated documentation by TravisCI"
git push

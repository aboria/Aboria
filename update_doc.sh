#!/bin/bash

ENCRYPTION_LABEL="b90db56f7413"
COMMIT_AUTHOR_EMAIL="martinjrobins@gmail.com"
        


SOURCE_BRANCH="master"
TARGET_BRANCH="gh-pages"

# Save some useful information
REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
SHA=`git rev-parse --verify HEAD`
HTML_DIR="doc/html-website"
BUILD_DIR=`pwd`

# Clone the existing gh-pages for this repo into out/
# Create a new empty branch if gh-pages doesn't exist yet (should only happen on first deply)
rm -rf $HTML_DIR
git clone $REPO $HTML_DIR
cd $HTML_DIR
git checkout $TARGET_BRANCH || git checkout --orphan $TARGET_BRANCH
cd $BUILD_DIR

# Clean out existing contents
rm -rf $HTML_DIR/**/* || exit 0

# Run our compile script
cmake .
make aboria-html-website

# Now let's go have some fun with the cloned repo
cd $HTML_DIR
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"

# If there are no changes to the compiled out (e.g. this is a README update) then just bail.
if git diff --quiet; then
    echo "No changes to the output on this push; exiting."
    exit 0
fi

# Commit the "changes", i.e. the new version.
# The delta will show diffs between new and old versions.
git add -A .
git commit -m "Deploy to GitHub Pages: ${SHA}"

# Get the deploy key by using Travis's stored variables to decrypt deploy_key.enc
ENCRYPTED_KEY_VAR="encrypted_${ENCRYPTION_LABEL}_key"
ENCRYPTED_IV_VAR="encrypted_${ENCRYPTION_LABEL}_iv"
ENCRYPTED_KEY=${!ENCRYPTED_KEY_VAR}
ENCRYPTED_IV=${!ENCRYPTED_IV_VAR}
openssl aes-256-cbc -K $ENCRYPTED_KEY -iv $ENCRYPTED_IV -in $BUILD_DIR/travis.enc -out travis -d
chmod 600 travis
eval `ssh-agent -s`
ssh-add travis

# Now that we're all set up, we can push.
git push $SSH_REPO $TARGET_BRANCH
eval `ssh-agent -k`


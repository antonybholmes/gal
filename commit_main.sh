# format code
#yarn update-version
#yarn format

# yarn run resize

# commit
#git add -A .
#git commit -m "Bug fixes and updates."
#git push -u origin dev

./commit.sh

git switch main
git merge dev
git push -u origin main

git switch dev
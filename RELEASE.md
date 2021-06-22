DENTIST Release Guide
=====================

![Document version: v2](https://img.shields.io/badge/Document%20version-v2-informational?logo=markdown)

This document describes the steps to take when preparing a new release.

## Quality Control

The following requirements must be met before creating a release:

1. The `develop` branch is checked out and the repo is in a clean state, i.e.
   `git status` is empty but for irrelevant untracked files.

2. The code compiles with `dub build`.

3. The pre-push hook `./test-current-commit.sh` succeeds. It includes a check
   for broken links, syntax check of the workflow and unit tests. Simply
   execute by pushing the `develop` branch.

4. Build Docker image with `docker build -t dentist:edge .`. This may take a
   while – go drink a coffee.

5. The example finishes successfully:
    
    ```sh
    # prepare example
    cd ../dentist-example
    # add the local repo if it is not yet a remote
    #git submodule foreach 'git remote add local "$(realpath "$PWD/../../../dentist")"'
    # update dentist submodule to latest develop
    git submodule foreach 'git fetch local && git checkout -f local/develop'
    # derive singularity image from Docker image
    make -B dentist_edge.sif
    # exhaustive tests using Singularity and pre-compiled binaries
    # go take a break (takes about 20min)
    make DENTIST_VERSION=edge test
    # reset dentist submodule
    git submodule update
    ```


## Prepare DENTIST release

The following steps must be taken to prepare the DENTIST release locally:

1. Update `CHANGELOG.md` (include correct version) and
   `git commit -m 'Updated CHANGELOG' CHANGELOG.md`.

2. Update `README.md`. Make sure to update all version numbers and release
   links. This can be changed later on but is best done before the release.

3. Bump version number with `./bump-version.sh`.

4. Make release tarball and Docker image with `./make-release.sh`.

5. Verify correct software version:
   ```sh
   DENTIST_VERSION=X.Y.Z
   docker run --rm aludi/dentist:"v$DENTIST_VERSION" dentist --version
   docker tag aludi/dentist:"v$DENTIST_VERSION" aludi/dentist:latest
   ```

6. Update `master` branch using `git checkout master && git merge --ff develop`


## Prepare example release

```sh
# prepare example
cd ../dentist-example

# update dentist submodule to the current release
git submodule foreach 'git fetch local && git checkout -f local/master'

# update DENTIST version in Makefile
sed -i -E 's/^(DENTIST_VERSION=).*/\1'"v$DENTIST_VERSION"'/' Makefile

# derive singularity image from Docker image
make "dentist_v$DENTIST_VERSION.sif"

# update DENTIST container in snakemake.yml
sed -i -E 's/^(dentist_container:).*/\1 '"dentist_v$DENTIST_VERSION.sif"'/' snakemake.yml

# exhaustive tests using Singularity and pre-compiled binaries
# go take a break (takes about 20min)
make test

# Update README; especially the version numbers
nano README.md

# commit changes
git add "dentist_v$DENTIST_VERSION.sif" bin Makefile README.md 

# clean up the worktree
rm dentist_edge.sif

# create the release tarball
make dist

# create release tag
git tag "v${DENTIST_VERSION}-1"
```


## Update GitHub Pages

```sh
# Switch to gh-pages
git checkout gh-pages
# Merge latest develop but stop before comitting
git merge --no-commit

# Resolve the merge...
git checkout gh-pages docs/_config.yml docs/_layouts/default.html \
  docs/index.md docs/Gemfile docs/Gemfile.lock
docker run --rm aludi/dentist:latest dentist --list-options > docs/list-of-commandline-options.md
# Prepare files for the webserver
( cd docs && ./prepare-site.sh )

# Create a combined commit
git merge --continue

# Go back to develop
git checkout develop
```


## Stage release

1. Publish DENTIST interally (staging):

    ```sh
    # Some links should be point to new release artifacts which are of course
    # not available yet. Ignore them for now.
    export IGNORE_LINKS="
        https://github.com/a-ludi/dentist/releases/download/v$DENTIST_VERSION/dentist.v$DENTIST_VERSION.x86_64.tar.gz
        https://github.com/a-ludi/dentist-example/releases/download/v$DENTIST_VERSION-1/dentist-example.tar.gz
    "

    # Push code changes to internal repository; this triggers the pre-push hook
    # which MUST succeed.
    git push internal develop master $(git tag --list 'v*')

    # Go ahead an check the internal repository website and see if everything
    # is in order.
    gio open https://git.mpi-cbg.de/ludwig/dentist

    # Create an internal release at
    gio open https://git.mpi-cbg.de/ludwig/dentist/-/releases/new

    # Make sure the release artifacts are present:
    [[ -f "dentist.v$DENTIST_VERSION.x86_64.tar.gz" ]] \
    || echo "error: missing release tarball"
    docker image ls aludi/dentist \
    | awk -v"version=v$DENTIST_VERSION" '
        ($2 == version) { iid = $3 }
        ($2 == "latest") { latest_iid = $3 }

        END {
            if (!iid)
                printf "error: release tag (%s) not found!\n", version;
            if (iid != latest_iid)
                printf "error: aludi/dentist:latest (%s) does not point to aludi/dentist:%s (%s). Please update using:\n    docker tag aludi/dentist:%s aludi/dentist:latest.\n", latest_iid, version, iid, version;

            exit !(iid && iid == latest_iid)
        }'
    ```

2. Publish the example interally (staging):

    ```sh
    cd ../dentist-example

    # Push code changes to internal repository
    git push internal master $(git tag --list 'v*')

    # Go ahead an check the internal repository website and see if everything
    # is in order. Check if the instructions are up-to-date and in sync with
    # DENTIST's README. Run the example with the online instructions (except
    # wget).
    gio open https://git.mpi-cbg.de/ludwig/dentist-example

    # Create an internal release at
    gio open https://git.mpi-cbg.de/ludwig/dentist-example/-/releases/new

    # Make sure the release artifacts are present:
    [[ -f "dentist_v$DENTIST_VERSION.sif" ]] \
    || echo "error: missing Singularity image"
    [[ -f "dentist-example.tar.gz" ]] \
    || echo "error: missing release tarball"
    [[ "dentist-example.tar.gz" -nt "../dentist_$DENTIST_VERSION.sif" ]] \
    || echo "error: release tarball is older than teh Singularity image"
    ```


## Publish release

**Attention:** This step makes the release public. Make sure to run all of the
above steps before publishing.

1. Publish DENTIST:

    ```sh
    # Push code changes to public repository; this triggers the pre-push hook
    # which MUST succeed.
    git push public develop master $(git tag --list 'v*')

    # Create a public release and remeber to include the release tarball
    # and the Singularity image. Use the title and text from the CHANGELOG.
    gio open https://github.com/a-ludi/dentist/releases/new

    # Publish Docker image
    docker push "aludi/dentist:v$DENTIST_VERSION"
    docker push aludi/dentist:latest

    # Update DUB package
    gio open https://code.dlang.org/my_packages/dentist
    ```

2. Publish DENTIST GitHub page:

    ```sh
    git push public --no-verify gh-pages
    ```

3. Publish the example:

    ```sh
    cd ../dentist-example

    # Push code changes to public repository
    git push public master $(git tag --list 'v*')

    # Create a public release and remeber to include the release tarball
    gio open https://github.com/a-ludi/dentist-example/releases/new
    ```

4. Push code changes to public repository:

    ```sh
    IGNORE_LINKS="
        https://github.com/a-ludi/dentist/releases/download/v$DENTIST_VERSION/dentist.v$DENTIST_VERSION.x86_64.tar.gz
        https://github.com/a-ludi/dentist-example/releases/download/v$DENTIST_VERSION-1/dentist-example.tar.gz" \
    git push public develop master $(git tag --list 'v*')
    ```


## Verify public release

1. Follow the instructions at <https://github.com/a-ludi/dentist/tree/develop#example>.

2. Check the GitHub page <https://a-ludi.github.io/dentist/>.

3. Check the DUB package at <https://code.dlang.org/packages/dentist>.

4. Check the Docker image at <https://hub.docker.com/repository/docker/aludi/dentist>.

DENTIST Release Guide
=====================

![Document version: v3](https://img.shields.io/badge/Document%20version-v3-informational?logo=markdown)

This document describes the steps to take when preparing a new release.

## Quality Control

The following requirements must be met before creating a release:

1. The `develop` branch is checked out and the repo is in a clean state, i.e.
   `git status` is empty but for irrelevant untracked files.

2. The code compiles with `dub build` and `dub build --config=testing`.

3. Unit tests pass with `dub test` and `dub test --config=testing`.

4. The pre-push hook `./test-current-commit.sh` succeeds. It includes a check
   for broken links, syntax check of the workflow, unit tests and integration
   tests. Simply execute by pushing the `develop` branch.

5. Build Conda packages locally:

    ```sh
    conda-build -c a_ludi -c bioconda \
        --clobber=conda/recipes/dentist-core/local-build.yaml \
        conda/recipes/dentist-core
    conda-build -c a_ludi -c bioconda -c conda-forge \
        --clobber=conda/recipes/dentist/local-build.yaml \
        conda/recipes/dentist
    ```

    This may take a while – go drink a coffee.

6. Build Singularity image:
    
    ```sh
    sudo singularity build example/dist/dentist_local.sif singularity/dentist_local.def
    ```

7. The example finishes successfully:
    
    ```sh
    # prepare example
    make -o dist/checksum.md5 -C example dentist_version=local prepare-dist

    # run using local Conda env
    make -o dist/checksum.md5 -C example dentist_version=local recreate-conda-env test-conda

    # run using pre-built binaries
    make -o dist/checksum.md5 -C example dentist_version=local test-binaries

    # run using Singularity image
    make -o dist/checksum.md5 -C example dentist_version=local test-singularity
    ```


## Prepare DENTIST release

The following steps must be taken to prepare the DENTIST release locally:

1. Update `CHANGELOG.md` (include correct version) and
   `git commit -m 'Updated CHANGELOG' CHANGELOG.md`.

2. Update `README.md` and `example/dist/README.md`. Make sure to update all
   version numbers and release links. Make sure the example section in
   `README.md` aggrees with `example/dist/README.md`. This can be changed later
   on but is best done before the release.

3. Update DENTIST version in `snakemake/Snakefile` (`dentist_container` and
   `dentist_env`). Create a new Singularity image definition and update the
   Conda recipes.

4. Check for obvious errors in the README: `./test-current-commit.sh`

5. Bump version number with `./bump-version.sh`.

6. Build Conda package locally:

    ```sh
    make dentist_version=local conda-dentist-core
    ```

7. Make release tarball and Docker image with `make dist`.

8. Verify correct software version:
   
    ```sh
    DENTIST_VERSION=vX.Y.Z
    dentist.$DENTIST_VERSION.x86_64/bin/dentist --version |& grep "$DENTIST_VERSION"
    ```

9. Update list of CLI options:

    ```sh
    make docs/list-of-commandline-options.md
    git commit -m 'Updated list of CLI options' docs/list-of-commandline-options.md
    ```

10. One more sanity check: `./test-current-commit.sh`

11. Update `master` branch using `git checkout master && git merge --ff develop`


## Update GitHub Pages

```sh
# Generate API docs
make -B docs.json

# Switch to gh-pages (second checkout recommended)
pushd ../dentist-gh-pages

# Copy API docs
cp ../dentist/docs.json ./docs.$DENTIST_VERSION.json

# Select new version as 'current'
ln -sf $DENTIST_VERSION api/current

# Generate HTML for API docs and update site files
make

# Commit changes
git add docs.$DENTIST_VERSION.json api/$DENTIST_VERSION
git commit -am "Updated site to DENTIST $DENTIST_VERSION"
git show --stat HEAD

# Go back to main repo
popd
```


## Stage release

Publish DENTIST internally (staging):

```sh
# Some links should point to new release artifacts which are of course
# not available yet. Ignore them for now.
export IGNORE_LINKS="
    https://github.com/a-ludi/dentist/releases/download/$DENTIST_VERSION/dentist.$DENTIST_VERSION.x86_64.tar.gz
    https://github.com/a-ludi/dentist/releases/download/$DENTIST_VERSION/dentist-example.tar.gz
"

# Push code changes to internal repository; this triggers the pre-push hook
# which MUST succeed.
git push internal master $(git tag --list 'v*')

# Go ahead an check the internal repository website and see if everything
# is in order.
gio open https://git.mpi-cbg.de/ludwig/dentist

# Create an internal release at
gio open https://git.mpi-cbg.de/ludwig/dentist/-/releases/new

# Make sure the release artifacts are present:
[[ -f "dentist.$DENTIST_VERSION.x86_64.tar.gz" ]] \
|| echo "error: missing release tarball"
```


## Publish release

**Attention:** This step makes the release public. Make sure to run all of the
above steps before publishing.

1. Publish DENTIST GitHub page:

    ```sh
    pushd ../dentist-gh-pages
    git push public gh-pages
    popd
    ```

2. Publish DENTIST:

    ```sh
    # Push code changes to public repository; this triggers the pre-push hook
    # which MUST succeed.
    git push public develop master $(git tag --list 'v*')

    # Create a public release draft and remember to include the release tarball
    # Use the title and text from the CHANGELOG.
    gio open https://github.com/a-ludi/dentist/releases/new

    # Update DUB package
    gio open https://code.dlang.org/my_packages/dentist
    ```

3. Create and publish Conda packages:

    ```sh
    # build and publish Conda packages
    pushd conda

    # Building the core package takes a while – go drink another coffee 
    ./conda-build.sh --upload recipes/dentist-core

    # consider building in another session as this takes a long time and is not
    # strictly required
    ./conda-build.sh --upload recipes/dentist

    # build also dependencies if they changed
    #./conda-build.sh --upload recipes/daccord
    #./conda-build.sh --upload recipes/daligner
    #./conda-build.sh --upload recipes/damapper
    #./conda-build.sh --upload recipes/damasker
    #./conda-build.sh --upload recipes/dascrubber
    #./conda-build.sh --upload recipes/dazz_db
    #./conda-build.sh --upload recipes/dentist
    #./conda-build.sh --upload recipes/dentist-core
    #./conda-build.sh --upload recipes/libunwind-static
    
    # back to repository root
    popd
    ```

4. Create and publish Singularity image:

    ```sh
    sudo singularity build dentist_${DENTIST_VERSION#v}.sif singularity/dentist_${DENTIST_VERSION#v}.def
    singularity sign dentist_${DENTIST_VERSION#v}.sif
    singularity push dentist_${DENTIST_VERSION#v}.sif library://a-ludi/default/dentist:${DENTIST_VERSION#v}
    ```

3. Create and publish the example:

    ```sh
    make -C example dist

    # Validate it works
    bash
        cd /tmp
        tar -xzf "$DENTIST/example/dentist-example.tar.gz"
        cd dentist-example/
        PATH="$PWD/bin:$PATH" snakemake --configfile=snakemake.yml --cores=all
        md5sum -c checksum.md5
        [[ -f "dentist_${DENTIST_VERSION#v}.sif" ]] || echo "error: SIF file missing"
    exit

    # Add example tarball to release and finish the draft release 
    gio open https://github.com/a-ludi/dentist-example/releases
    ```

## Verify public release

1. Follow the instructions at <https://github.com/a-ludi/dentist#example>.

2. Check the GitHub page <https://a-ludi.github.io/dentist/>.

3. Check the DUB package at <https://code.dlang.org/packages/dentist>.

4. Check the Docker image at <https://hub.docker.com/repository/docker/aludi/dentist>.

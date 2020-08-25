#!/bin/sh

echo "refs/dummy" | "$(dirname "$0")/.githooks/pre-push" localhost https://localhost

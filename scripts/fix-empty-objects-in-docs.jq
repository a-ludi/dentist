. as $docs
| paths(
    objects
    | select(has("members"))
    | .members
    | map(select(has("name") | not))
    | select(length > 0)
)
| . as $path
| $docs
| getpath($path)
| .members
| map(select(has("name")))
| . as $fixed
| $docs
| setpath([$path[], "members"]; $fixed)

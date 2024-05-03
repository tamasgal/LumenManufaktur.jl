# An example

The following function determines the meaning of life.

```@example usage
using LumenManufaktur

meaningoflife()
```

Examples with the same "tag" (like `usage` above) share the same Julia
process, so that everything is in the same scope. The package is therefore already
imported, so we can determine the meaning of life again `;)`

```@example usage
meaningoflife()
```

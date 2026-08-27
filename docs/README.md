# Building these docs

Moved here from the documentation home page, and corrected: with
DocumenterVitepress the build no longer produces a browsable
`docs/build/index.html`.

From the repository root:

```
julia --project=docs docs/make.jl
```

The site lands in `docs/build/1`, and its links have no `.html` suffix
(`cleanUrls`), so a plain static server answers 404 to all of them. Serve it
with VitePress instead, from `docs/`:

```
npm run docs:preview
```

`npm` comes from the `NodeJS_20_jll` artifact that DocumenterVitepress pulls in,
so nothing has to be installed system-wide; if it is not on `PATH`:

```
PATH="$(dirname $(ls ~/.julia/artifacts/*/bin/npm | head -1)):$PATH" npm run docs:preview
```

`npm run docs:dev` does the same with hot reload while editing.

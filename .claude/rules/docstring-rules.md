---
paths:
  - src/**/*.jl
---

# Docstring Rules

## Use DocStringExtensions.jl

Pick the abbreviation by what the docstring is attached to:

- **Functions / methods / constructors:** Always use `$(TYPEDSIGNATURES)` for the 
  automatic type-annotated signature.
- **Structs / types:** `$(TYPEDEF)` for the type definition line, and `$(TYPEDFIELDS)`
  for the field list.
- Document struct fields with a **string docstring above each field** as `$(TYPEDFIELDS)`
  only picks up string docstrings (the inline `#` comments are ignored).

Two cases where `$(TYPEDSIGNATURES)` does not work — keep a **manual signature line** instead:

- **Interface stubs** documented on a method-less `function foo end`.
- **One docstring covering several call forms** (e.g. a constructor with two signatures).
  `$(TYPEDSIGNATURES)` shows only the one method it is attached to; a manual list keeps all
  forms visible.

## CRITICAL: Always use `jldoctest`, NEVER plain `julia` blocks

Plain code blocks (`` ```julia ``) are NOT tested and can become stale or incorrect.
Doctests (`` ```jldoctest ``) are automatically tested and verified to work.

### Example:

~~~~
"""
$(TYPEDEF)

Grid definition with a given floating point precision, size and snow/soil layer
thicknesses. Construct with `Grid(Tf; kwargs...)`:

```jldoctest
using FlexibleSnowModelOSHD

grid = Grid(Float32, Nx = 10, Ny = 5)

# output
Grid
├── Precision: Float32
├── Nx: 10
├── Ny: 5
├── Dzsnow: [0.1, 0.2, 0.4]
└── Dzsoil: [0.1, 0.2, 0.4, 0.8]
```

# Fields

$(TYPEDFIELDS)
"""
~~~~

## Doctest Best Practices

- Always include expected output after `# output`
- Use simple, verifiable output (e.g., `typeof(result)`, accessing a field)
- Doctests should exercise `Base.show` to verify objects display correctly
- Keep doctests minimal but complete enough to verify the feature works
- **Do NOT use boolean comparisons as the final line** (e.g., avoid `x ≈ 1.0` or `obj isa Type`)
- Instead, make the final line invoke a `show` method that prints something useful

# io/FloodA5Keys.jl — access keys and tokens (Cesium ion, Google Map Tiles, …)
#
# A key is looked up, in order, in
#   1. the environment variable of that name (CI, servers, containers);
#   2. the file named by the environment variable FLOODA5_KEYS, if set;
#   3. keys.json in the repository root (git-ignored);
#   4. ~/.flooda5/keys.json (per user: one copy for every clone and worktree).
# The files are JSON objects keyed by the same names (template:
# keys.example.json in the repository root). Empty values and template
# placeholders ("YOUR_…") count as unset. Values are never logged:
# `key_source` says where a key came from without showing it.

module FloodA5Keys

using JSON3

export flooda5_key, key_source

"The keys FloodA5 knows (keys.example.json lists the same)."
const KEY_NAMES = ("FLOODA5_CESIUM_ION_TOKEN", "FLOODA5_GOOGLE_TILES_KEY")

repo_keys_file() = joinpath(dirname(@__DIR__), "keys.json")
user_keys_file() = joinpath(homedir(), ".flooda5", "keys.json")

"The key files searched, in order (FLOODA5_KEYS, the repository's, the user's)."
function keys_files()
    fs = String[]
    k = strip(get(ENV, "FLOODA5_KEYS", ""))
    isempty(k) || push!(fs, String(k))
    push!(fs, repo_keys_file(), user_keys_file())
    return fs
end

_usable(v) = v isa AbstractString && !isempty(strip(v)) && !startswith(strip(v), "YOUR_")

const _WARNED = Set{String}()

# The file's keys, or an empty Dict (with one warning, naming the file only).
function _read_keys(f::AbstractString)
    isfile(f) || return Dict{String,Any}()
    try
        # (a UTF-8 byte-order mark, as Notepad and Windows PowerShell 5.1 write, is dropped)
        return JSON3.read(lstrip(==('﻿'), read(f, String)), Dict{String,Any})
    catch
        f in _WARNED || (push!(_WARNED, f); @warn "FloodA5Keys: $f is not a JSON object of keys; ignoring it")
        return Dict{String,Any}()
    end
end

"""
    key_lookup(name; files = keys_files()) → (value, source) or (nothing, nothing)

`source` is "environment" or the path of the file the key came from.
"""
function key_lookup(name::AbstractString; files = keys_files())
    v = get(ENV, name, "")
    _usable(v) && return (String(strip(v)), "environment")
    for f in files
        v = get(_read_keys(f), name, nothing)
        _usable(v) && return (String(strip(v)), f)
    end
    return (nothing, nothing)
end

"""
    flooda5_key(name; default = "", files = keys_files()) → String

The key `name` from the environment or the key files (see the module comment),
or `default`.
"""
flooda5_key(name::AbstractString; default = "", files = keys_files()) =
    something(key_lookup(name; files)[1], default)

"""
    key_source(name; files = keys_files()) → String or nothing

Where `name` comes from ("environment" or a file path), without its value.
"""
key_source(name::AbstractString; files = keys_files()) = key_lookup(name; files)[2]

end # module FloodA5Keys

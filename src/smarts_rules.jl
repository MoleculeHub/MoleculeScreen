# Load SMARTS rules from CSV
function load_smarts_rules()
    csv_path = joinpath(@__DIR__, "..", "data", "rules_data.csv")

    if !isfile(csv_path)
        @warn "SMARTS rules CSV file not found at $csv_path. Using empty rule set."
        return (Dict{String, Vector{String}}(), Dict{String, Vector{String}}())
    end

    try
        df = CSV.read(csv_path, DataFrame)

        rules = Dict{String, Vector{String}}()
        descriptions = Dict{String, Vector{String}}()

        for row in eachrow(df)
            rule_set = lowercase(String(row.rule_set_name))
            smarts = String(row.smarts)
            description = String(row.description)

            if haskey(rules, rule_set)
                push!(rules[rule_set], smarts)
                push!(descriptions[rule_set], description)
            else
                rules[rule_set] = [smarts]
                descriptions[rule_set] = [description]
            end
        end

        return (rules, descriptions)
    catch e
        @warn "Error loading SMARTS rules: $e. Using empty rule set."
        return (Dict{String, Vector{String}}(), Dict{String, Vector{String}}())
    end
end

const SMARTS_RULES, SMARTS_DESCRIPTIONS = load_smarts_rules()

"""
    get_available_smarts_filters() -> Vector{String}

Get a list of all available SMARTS filter names loaded from the CSV file.

# Returns

  - `Vector{String}`: List of available filter names
"""
function get_available_smarts_filters()
    return collect(keys(SMARTS_RULES))
end


"""
    check_smarts_filter(mol::Molecule, filter_name::String) -> Bool

Check if a molecule passes a specific SMARTS-based filter.

# Arguments

  - `mol::Molecule`: Molecule to check
  - `filter_name::String`: Name of the filter to apply (case-insensitive)

Available filters (loaded from data/rules_data.csv):

  - "pains" - Pan-Assay Interference Compounds
  - "elililly" - Eli Lilly filters
  - "surechembl" - SureChEMBL reactivity alerts
  - "mlsmr" - MLSMR filters
  - "bms" - Bristol-Myers Squibb filters
  - "brenk" - Brenk filter for unwanted functional groups
  - "inpharmatica" - Inpharmatica filters
  - "glaxo" - GlaxoSmithKline filters
  - "lint" - Lilly LINT filters

Use `get_available_smarts_filters()` to see all loaded filter names.

# Returns

  - `Bool`: true if molecule passes (no alerts), false if molecule contains unwanted substructures
"""
function check_smarts_filter(mol::Molecule, filter_name::String)
    !mol.valid && return false

    # for case-insensitive matching
    filter_key = lowercase(filter_name)

    if !haskey(SMARTS_RULES, filter_key)
        available_filters = join(sort(collect(keys(SMARTS_RULES))), ", ")
        @warn "Unknown SMARTS filter: '$filter_name'. Available filters: $available_filters"
        return false
    end

    patterns = SMARTS_RULES[filter_key]

    for pattern in patterns
        try
            if has_substructure_match(mol, pattern)
                return false
            end
        catch e
            @warn "Error checking pattern '$pattern' in filter '$filter_name': $e"
            continue
        end
    end

    return true
end

"""
    apply_smarts_filters(mol::Molecule; filters::Vector{String}=["pains"]) -> Dict{String, Any}

Apply multiple SMARTS-based filters to a molecule.

# Arguments

  - `mol::Molecule`: Molecule to filter
  - `filters::Vector{String}`: List of SMARTS filters to apply

# Returns

  - `Dict{String, Any}`: Dictionary mapping filter names to results with pass/fail status and violation descriptions
"""
function apply_smarts_filters(mol::Molecule; filters::Vector{String} = ["pains"])
    results = Dict{String, Any}()

    for filter_name in filters
        pass_status = check_smarts_filter(mol, filter_name)
        violations = get_smarts_violations_with_descriptions(mol, filter_name)

        results[filter_name] = Dict(
            "pass" => pass_status,
            "violations" => violations
        )
    end

    return results
end

"""
    get_smarts_violations(mol::Molecule, filter_name::String) -> Vector{String}

Get the specific SMARTS patterns that match in a molecule for a given filter.

# Arguments

  - `mol::Molecule`: Molecule to check
  - `filter_name::String`: Name of the filter to check

# Returns

  - `Vector{String}`: List of SMARTS patterns that match (violations)
"""
function get_smarts_violations(mol::Molecule, filter_name::String)
    !mol.valid && return String[]

    filter_key = lowercase(filter_name)

    if !haskey(SMARTS_RULES, filter_key)
        available_filters = join(sort(collect(keys(SMARTS_RULES))), ", ")
        @warn "Unknown SMARTS filter: '$filter_name'. Available filters: $available_filters"
        return String[]
    end

    patterns = SMARTS_RULES[filter_key]
    violations = String[]

    for pattern in patterns
        try
            if has_substructure_match(mol, pattern)
                push!(violations, pattern)
            end
        catch e
            @warn "Error checking pattern '$pattern' in filter '$filter_name': $e"
            continue
        end
    end

    return violations
end

"""
    get_smarts_violations_with_descriptions(mol::Molecule, filter_name::String) -> Vector{NamedTuple{(:smarts, :description), Tuple{String, String}}}

Get the specific SMARTS patterns that match in a molecule along with their descriptions.

# Arguments

  - `mol::Molecule`: Molecule to check
  - `filter_name::String`: Name of the filter to check

# Returns

  - `Vector{NamedTuple}`: List of (smarts=pattern, description=desc) named tuples for violations
"""
function get_smarts_violations_with_descriptions(mol::Molecule, filter_name::String)
    !mol.valid && return NamedTuple{(:smarts, :description), Tuple{String, String}}[]

    filter_key = lowercase(filter_name)

    if !haskey(SMARTS_RULES, filter_key)
        available_filters = join(sort(collect(keys(SMARTS_RULES))), ", ")
        @warn "Unknown SMARTS filter: '$filter_name'. Available filters: $available_filters"
        return NamedTuple{(:smarts, :description), Tuple{String, String}}[]
    end

    patterns = SMARTS_RULES[filter_key]
    descriptions = SMARTS_DESCRIPTIONS[filter_key]
    violations = NamedTuple{(:smarts, :description), Tuple{String, String}}[]

    for (pattern, description) in zip(patterns, descriptions)
        try
            if has_substructure_match(mol, pattern)
                push!(violations, (smarts=pattern, description=description))
            end
        catch e
            @warn "Error checking pattern '$pattern' in filter '$filter_name': $e"
            continue
        end
    end

    return violations
end

"""
    is_clean_molecule(mol::Molecule; smarts_filters::Vector{String}=["pains", "brenk"]) -> Bool

Check if a molecule is "clean" (passes all specified SMARTS filters).

# Arguments

  - `mol::Molecule`: Molecule to check
  - `smarts_filters::Vector{String}`: List of SMARTS filters to apply

# Returns

  - `Bool`: true if molecule passes all filters, false otherwise
"""
function is_clean_molecule(
    mol::Molecule; smarts_filters::Vector{String} = ["pains", "brenk"]
)
    results = apply_smarts_filters(mol; filters = smarts_filters)
    return all(r["pass"] for r in values(results))
end

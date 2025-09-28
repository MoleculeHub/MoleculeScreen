using MoleculeFlow

# Property-based filtering functions for drug-like compounds

"""
    lipinski_ro5(mol::Molecule) -> Bool

Apply Lipinski's Rule of Five filter.

Compounds should satisfy at least 3 of 4 criteria:
- Molecular weight ≤ 500 Da
- LogP ≤ 5
- H-bond donors ≤ 5
- H-bond acceptors ≤ 10

# Arguments
- `mol::Molecule`: Molecule to filter

# Returns
- `Bool`: true if molecule passes (≤ 1 violation), false otherwise
"""
function lipinski_ro5(mol::Molecule)
    !mol.valid && return false

    violations = 0

    # MW ≤ 500 Da
    violations += molecular_weight(mol) > 500 ? 1 : 0

    # logP ≤ 5
    violations += logp(mol) > 5 ? 1 : 0

    # H-bond donors ≤ 5
    violations += num_hbd(mol) > 5 ? 1 : 0

    # H-bond acceptors ≤ 10
    violations += num_hba(mol) > 10 ? 1 : 0

    return violations <= 1
end

"""
    veber_rules(mol::Molecule) -> Bool

Apply Veber's rules for oral bioavailability.

Criteria:
- Rotatable bonds ≤ 10
- Topological polar surface area (TPSA) ≤ 140 Å²

# Arguments
- `mol::Molecule`: Molecule to filter

# Returns
- `Bool`: true if molecule passes both criteria, false otherwise
"""
function veber_rules(mol::Molecule)
    !mol.valid && return false

    return num_rotatable_bonds(mol) <= 10 && tpsa(mol) <= 140
end

"""
    ghose_filter(mol::Molecule) -> Bool

Apply Ghose filter for drug-like chemical space.

Criteria:
- MW: 160–480 Da
- LogP: -0.4 to +5.6
- Molar refractivity: 40–130
- Heavy atom count: 20–70

# Arguments
- `mol::Molecule`: Molecule to filter

# Returns
- `Bool`: true if molecule passes all criteria, false otherwise
"""
function ghose_filter(mol::Molecule)
    !mol.valid && return false

    mw = molecular_weight(mol)
    lp = logp(mol)
    mr = molecular_weight(mol) * 0.3
    hac = heavy_atom_count(mol)

    return (160 <= mw <= 480) && (-0.4 <= lp <= 5.6) && (40 <= mr <= 130) && (20 <= hac <= 70)
end

"""
    egan_filter(mol::Molecule) -> Bool

Apply Egan's filter for passive permeability.

Criteria:
- LogP ≤ 5.88
- TPSA ≤ 131 Å²

# Arguments
- `mol::Molecule`: Molecule to filter

# Returns
- `Bool`: true if molecule passes both criteria, false otherwise
"""
function egan_filter(mol::Molecule)
    !mol.valid && return false

    return logp(mol) <= 5.88 && tpsa(mol) <= 131
end

"""
    muegge_filter(mol::Molecule) -> Bool

Apply Muegge filter for general drug-likeness.

Criteria:
- MW: 200–600 Da
- LogP: -2 to 5
- TPSA ≤ 150 Å²
- Rotatable bonds ≤ 15
- H-bond donors ≤ 5
- H-bond acceptors ≤ 10

# Arguments
- `mol::Molecule`: Molecule to filter

# Returns
- `Bool`: true if molecule passes all criteria, false otherwise
"""
function muegge_filter(mol::Molecule)
    !mol.valid && return false

    mw = molecular_weight(mol)
    lp = logp(mol)
    tps = tpsa(mol)
    rb = num_rotatable_bonds(mol)
    hbd = num_hbd(mol)
    hba = num_hba(mol)

    return (200 <= mw <= 600) && (-2 <= lp <= 5) && (tps <= 150) && (rb <= 15) && (hbd <= 5) && (hba <= 10)
end

"""
    pfizer_3_75_rule(mol::Molecule) -> Bool

Apply Pfizer's 3/75 rule to identify compounds prone to safety issues.

Compounds are flagged as potentially problematic if:
- LogP > 3 AND TPSA < 75 Å²

# Arguments
- `mol::Molecule`: Molecule to filter

# Returns
- `Bool`: true if molecule passes (not flagged), false if flagged as problematic
"""
function pfizer_3_75_rule(mol::Molecule)
    !mol.valid && return false

    # Return false if compound is flagged (LogP > 3 AND TPSA < 75)
    return !(logp(mol) > 3 && tpsa(mol) < 75)
end

"""
    gsk_4_400_rule(mol::Molecule) -> Bool

Apply GSK's 4/400 rule for CNS drug quality.

Criteria:
- H-bond donors ≤ 4
- MW ≤ 400 Da

# Arguments
- `mol::Molecule`: Molecule to filter

# Returns
- `Bool`: true if molecule passes both criteria, false otherwise
"""
function gsk_4_400_rule(mol::Molecule)
    !mol.valid && return false

    return num_hbd(mol) <= 4 && molecular_weight(mol) <= 400
end

"""
    golden_triangle(mol::Molecule) -> Bool

Apply Golden Triangle filter for optimal ADMET properties.

Criteria:
- LogP: -0.5 to 4.5
- MW: 200–500 Da

# Arguments
- `mol::Molecule`: Molecule to filter

# Returns
- `Bool`: true if molecule is in the golden triangle, false otherwise
"""
function golden_triangle(mol::Molecule)
    !mol.valid && return false

    lp = logp(mol)
    mw = molecular_weight(mol)

    return (-0.5 <= lp <= 4.5) && (200 <= mw <= 500)
end

"""
    apply_property_filters(mol::Molecule; filters::Vector{Symbol}=[:lipinski]) -> Dict{Symbol, Bool}

Apply multiple property-based filters to a molecule.

# Arguments
- `mol::Molecule`: Molecule to filter
- `filters::Vector{Symbol}`: List of filters to apply

Available filters:
- `:lipinski` - Lipinski's Rule of Five
- `:veber` - Veber's rules
- `:ghose` - Ghose filter
- `:egan` - Egan's filter
- `:muegge` - Muegge filter
- `:pfizer` - Pfizer's 3/75 rule
- `:gsk` - GSK's 4/400 rule
- `:golden_triangle` - Golden Triangle

# Returns
- `Dict{Symbol, Bool}`: Dictionary mapping filter names to results
"""
function apply_property_filters(mol::Molecule; filters::Vector{Symbol}=[:lipinski])
    filter_functions = Dict(
        :lipinski => lipinski_ro5,
        :veber => veber_rules,
        :ghose => ghose_filter,
        :egan => egan_filter,
        :muegge => muegge_filter,
        :pfizer => pfizer_3_75_rule,
        :gsk => gsk_4_400_rule,
        :golden_triangle => golden_triangle
    )

    results = Dict{Symbol, Bool}()
    for filter_name in filters
        if haskey(filter_functions, filter_name)
            results[filter_name] = filter_functions[filter_name](mol)
        else
            @warn "Unknown filter: $filter_name"
        end
    end

    return results
end
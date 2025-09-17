@testset "Filtering Tests" begin

    # Setup test data
    function setup_filtering_test_data()
        alerts = [
            StructuralAlert(1, 1, "Alcohol", "[OH]", "Basic", 8, 0),
            StructuralAlert(2, 1, "Ketone", "C=O", "Basic", 7, 0),
            StructuralAlert(3, 1, "Benzene", "c1ccccc1", "Basic", 6, 0),
            StructuralAlert(4, 2, "Amine", "[NH2]", "Nitrogen", 5, 0),
            StructuralAlert(5, 2, "Pyridine", "c1ccncc1", "Nitrogen", 8, 0)
        ]
        collection = AlertCollection(alerts)

        molecules = [
            mol_from_smiles("CCO"),        # Ethanol - alcohol (high priority)
            mol_from_smiles("CC(=O)C"),    # Acetone - ketone (high priority)
            mol_from_smiles("c1ccccc1"),   # Benzene - benzene (low priority)
            mol_from_smiles("CCN"),        # Ethylamine - amine (low priority)
            mol_from_smiles("CCC"),        # Propane - no alerts
            mol_from_smiles("c1ccncc1")    # Pyridine - pyridine (high priority)
        ]

        mol_ids = ["ethanol", "acetone", "benzene", "ethylamine", "propane", "pyridine"]

        return collection, molecules, mol_ids
    end

    @testset "filter_molecules_by_alerts - Return Passed" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        passed_mols, passed_ids, results = filter_molecules_by_alerts(
            molecules, collection, molecule_ids=mol_ids, priority_threshold=7, return_passed=true
        )

        # Should return molecules with no high-priority alerts (priority < 7)
        # Benzene (priority 6), Ethylamine (priority 5), Propane (no alerts)
        @test length(passed_mols) == 3
        @test "benzene" in passed_ids
        @test "ethylamine" in passed_ids
        @test "propane" in passed_ids
        @test length(results) == 3
        @test all(r.passed for r in results)
    end

    @testset "filter_molecules_by_alerts - Return Failed" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        failed_mols, failed_ids, results = filter_molecules_by_alerts(
            molecules, collection, molecule_ids=mol_ids, priority_threshold=7, return_passed=false
        )

        # Should return molecules with high-priority alerts (priority >= 7)
        # Ethanol (priority 8), Acetone (priority 7), Pyridine (priority 8)
        @test length(failed_mols) == 3
        @test "ethanol" in failed_ids
        @test "acetone" in failed_ids
        @test "pyridine" in failed_ids
        @test length(results) == 3
        @test all(!r.passed for r in results)
    end

    @testset "filter_molecules_by_alerts - SMILES Input" begin
        collection, _, mol_ids = setup_filtering_test_data()
        smiles_list = ["CCO", "CC(=O)C", "CCC"]  # Ethanol, Acetone, Propane

        passed_smiles, passed_ids, results = filter_molecules_by_alerts(
            smiles_list, collection, molecule_ids=mol_ids[1:3], priority_threshold=7, return_passed=true
        )

        # Only Propane should pass (no alerts)
        @test length(passed_smiles) == 1
        @test passed_smiles[1] == "CCC"
        @test passed_ids[1] == "propane"
    end

    @testset "filter_molecules_by_alerts - Auto-generated IDs" begin
        collection, molecules, _ = setup_filtering_test_data()

        passed_mols, passed_ids, results = filter_molecules_by_alerts(
            molecules[1:3], collection, priority_threshold=7, return_passed=true
        )

        # Should have auto-generated IDs
        @test all(startswith(id, "mol_") for id in passed_ids)
        @test length(passed_mols) == 1  # Only benzene passes
    end

    @testset "filter_by_specific_alerts" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        # Filter out molecules with alcohol and ketone alerts (rule IDs 1 and 2)
        clean_mols, clean_ids, results = filter_by_specific_alerts(
            molecules, collection, [1, 2], molecule_ids=mol_ids
        )

        # Should exclude ethanol (has alcohol) and acetone (has ketone)
        @test length(clean_mols) == 4
        @test "ethanol" ∉ clean_ids
        @test "acetone" ∉ clean_ids
        @test "benzene" in clean_ids
        @test "ethylamine" in clean_ids
        @test "propane" in clean_ids
        @test "pyridine" in clean_ids

        # Check results
        @test length(results) == 6  # All molecules have results
        ethanol_result = results[findfirst(r -> r.molecule_id == "ethanol", results)]
        @test !ethanol_result.passed
        @test 1 in ethanol_result.alert_matches

        propane_result = results[findfirst(r -> r.molecule_id == "propane", results)]
        @test propane_result.passed
        @test isempty(propane_result.alert_matches)
    end

    @testset "filter_by_specific_alerts - Nonexistent Alert ID" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        # Try to filter by nonexistent alert ID
        clean_mols, clean_ids, results = filter_by_specific_alerts(
            molecules, collection, [999], molecule_ids=mol_ids
        )

        # All molecules should pass since alert doesn't exist
        @test length(clean_mols) == 6
        @test length(clean_ids) == 6
    end

    @testset "filter_by_rule_set" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        # Filter using only "Basic" rule set
        basic_filtered, basic_ids, results = filter_by_rule_set(
            molecules, collection, "Basic", molecule_ids=mol_ids, priority_threshold=7
        )

        # Basic rule set has: Alcohol (8), Ketone (7), Benzene (6)
        # So ethanol and acetone should fail, others should pass
        @test length(basic_filtered) == 4
        @test "ethanol" ∉ basic_ids  # Has alcohol (priority 8)
        @test "acetone" ∉ basic_ids  # Has ketone (priority 7)
        @test "benzene" in basic_ids     # Benzene (priority 6) passes
        @test "ethylamine" in basic_ids  # No Basic alerts
        @test "propane" in basic_ids     # No alerts
        @test "pyridine" in basic_ids    # No Basic alerts
    end

    @testset "filter_by_rule_set - Nonexistent Rule Set" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        filtered_mols, filtered_ids, results = filter_by_rule_set(
            molecules, collection, "NonExistent", molecule_ids=mol_ids
        )

        # All molecules should pass since no alerts from this rule set
        @test length(filtered_mols) == 6
        @test all(r.passed for r in results)
    end

    @testset "filter_by_functional_groups - Exclude Groups" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        # Exclude molecules with alcohol functional group
        clean_mols, clean_ids = filter_by_functional_groups(
            molecules, [:alcohol], molecule_ids=mol_ids, exclude_groups=true
        )

        # Should exclude ethanol (has alcohol)
        @test "ethanol" ∉ clean_ids
        @test "acetone" in clean_ids
        @test "propane" in clean_ids
        @test length(clean_mols) == 5
    end

    @testset "filter_by_functional_groups - Keep Only Groups" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        # Keep only molecules with alcohol functional group
        alcohol_mols, alcohol_ids = filter_by_functional_groups(
            molecules, [:alcohol], molecule_ids=mol_ids, exclude_groups=false
        )

        # Should keep only ethanol
        @test length(alcohol_mols) == 1
        @test alcohol_ids[1] == "ethanol"
    end

    @testset "filter_by_functional_groups - Multiple Groups" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        # Exclude molecules with alcohol or ketone groups
        clean_mols, clean_ids = filter_by_functional_groups(
            molecules, [:alcohol, :ketone], molecule_ids=mol_ids, exclude_groups=true
        )

        # Should exclude ethanol (alcohol) and acetone (ketone)
        @test "ethanol" ∉ clean_ids
        @test "acetone" ∉ clean_ids
        @test length(clean_mols) == 4
    end

    @testset "filter_by_functional_groups - Auto-generated IDs" begin
        collection, molecules, _ = setup_filtering_test_data()

        clean_mols, clean_ids = filter_by_functional_groups(
            molecules[1:3], [:alcohol], exclude_groups=true
        )

        @test all(startswith(id, "mol_") for id in clean_ids)
        @test length(clean_mols) == 2  # Exclude ethanol
    end

    @testset "create_custom_filter" begin
        collection, _, _ = setup_filtering_test_data()

        # Create custom filter with multiple criteria
        criteria = Dict{String, Any}(
            "rule_sets" => ["Basic"],
            "min_priority" => 7,
            "exclude_rule_ids" => [3]  # Exclude benzene
        )

        custom_collection = create_custom_filter(collection, criteria)

        # Should have alcohol and ketone from Basic rule set with priority >= 7, excluding benzene
        @test length(custom_collection.alerts) == 2
        alert_ids = [alert.rule_id for alert in custom_collection.alerts]
        @test 1 in alert_ids  # Alcohol
        @test 2 in alert_ids  # Ketone
        @test 3 ∉ alert_ids   # Benzene excluded
        @test 4 ∉ alert_ids   # Amine not in Basic rule set
        @test 5 ∉ alert_ids   # Pyridine not in Basic rule set
    end

    @testset "create_custom_filter - Rule Sets Only" begin
        collection, _, _ = setup_filtering_test_data()

        criteria = Dict{String, Any}("rule_sets" => ["Nitrogen"])
        custom_collection = create_custom_filter(collection, criteria)

        @test length(custom_collection.alerts) == 2
        rule_set_names = [alert.rule_set_name for alert in custom_collection.alerts]
        @test all(name == "Nitrogen" for name in rule_set_names)
    end

    @testset "create_custom_filter - Priority Range" begin
        collection, _, _ = setup_filtering_test_data()

        criteria = Dict{String, Any}(
            "min_priority" => 6,
            "max_priority" => 7
        )
        custom_collection = create_custom_filter(collection, criteria)

        @test length(custom_collection.alerts) == 2  # Ketone (7) and Benzene (6)
        priorities = [alert.priority for alert in custom_collection.alerts]
        @test all(6 <= p <= 7 for p in priorities)
    end

    @testset "create_custom_filter - Specific Rule IDs" begin
        collection, _, _ = setup_filtering_test_data()

        criteria = Dict{String, Any}("rule_ids" => [1, 3, 5])
        custom_collection = create_custom_filter(collection, criteria)

        @test length(custom_collection.alerts) == 3
        alert_ids = [alert.rule_id for alert in custom_collection.alerts]
        @test alert_ids == [1, 3, 5]
    end

    @testset "create_custom_filter - Empty Result" begin
        collection, _, _ = setup_filtering_test_data()

        # Impossible criteria
        criteria = Dict{String, Any}(
            "rule_sets" => ["NonExistent"],
            "min_priority" => 10
        )
        custom_collection = create_custom_filter(collection, criteria)

        @test length(custom_collection.alerts) == 0
    end

    @testset "Filtering Edge Cases" begin
        collection, molecules, mol_ids = setup_filtering_test_data()

        # Test with empty molecule list
        empty_mols, empty_ids, empty_results = filter_molecules_by_alerts(
            Molecule[], collection, molecule_ids=String[]
        )
        @test length(empty_mols) == 0
        @test length(empty_ids) == 0
        @test length(empty_results) == 0

        # Test with invalid molecules
        invalid_mol = Molecule(_rdkit_mol=nothing, valid=false, source="invalid")
        invalid_filtered, invalid_ids, invalid_results = filter_molecules_by_alerts(
            [invalid_mol], collection, molecule_ids=["invalid"]
        )
        @test length(invalid_filtered) == 0  # Invalid molecules don't pass
        @test length(invalid_results) == 1
        @test !invalid_results[1].passed
    end
end
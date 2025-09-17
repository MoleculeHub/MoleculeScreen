@testset "Screening Tests" begin

    # Setup test molecules and alerts
    function setup_test_data()
        alerts = [
            StructuralAlert(1, 1, "Alcohol", "[OH]", "Basic", 8, 0),
            StructuralAlert(2, 1, "Ketone", "C=O", "Basic", 7, 0),
            StructuralAlert(3, 1, "Benzene", "c1ccccc1", "Basic", 6, 0),
            StructuralAlert(4, 1, "Amine", "[NH2]", "Basic", 5, 0)
        ]
        collection = AlertCollection(alerts)

        molecules = [
            mol_from_smiles("CCO"),        # Ethanol - has alcohol
            mol_from_smiles("CC(=O)C"),    # Acetone - has ketone
            mol_from_smiles("c1ccccc1"),   # Benzene - has benzene
            mol_from_smiles("CCN"),        # Ethylamine - has amine
            mol_from_smiles("CCC"),        # Propane - no alerts
            mol_from_smiles("CC(=O)c1ccccc1") # Acetophenone - has ketone and benzene
        ]

        return collection, molecules
    end

    @testset "ScreeningResult Construction" begin
        result = ScreeningResult("test_mol", 3, [1, 2, 5], 2, false)
        @test result.molecule_id == "test_mol"
        @test result.total_alerts == 3
        @test result.alert_matches == [1, 2, 5]
        @test result.high_priority_count == 2
        @test result.passed == false
    end

    @testset "screen_molecule - Single Molecule" begin
        collection, molecules = setup_test_data()

        # Test ethanol (should match alcohol alert)
        result = screen_molecule(molecules[1], collection, molecule_id="ethanol")
        @test result.molecule_id == "ethanol"
        @test result.total_alerts == 1
        @test 1 in result.alert_matches  # Alcohol alert
        @test result.high_priority_count == 1  # Alcohol has priority 8
        @test result.passed == false  # Has high priority alert

        # Test propane (should have no alerts)
        result = screen_molecule(molecules[5], collection, molecule_id="propane")
        @test result.molecule_id == "propane"
        @test result.total_alerts == 0
        @test isempty(result.alert_matches)
        @test result.high_priority_count == 0
        @test result.passed == true  # No alerts

        # Test acetophenone (should match ketone and benzene)
        result = screen_molecule(molecules[6], collection, molecule_id="acetophenone")
        @test result.total_alerts == 2
        @test 2 in result.alert_matches  # Ketone
        @test 3 in result.alert_matches  # Benzene
        @test result.high_priority_count == 1  # Only ketone has priority >= 7
    end

    @testset "screen_molecule - SMILES String Input" begin
        collection, _ = setup_test_data()

        result = screen_molecule("CCO", collection, molecule_id="ethanol_smiles")
        @test result.molecule_id == "ethanol_smiles"
        @test result.total_alerts == 1
        @test 1 in result.alert_matches  # Alcohol alert
    end

    @testset "screen_molecule - Invalid Molecule" begin
        collection, _ = setup_test_data()

        # Create invalid molecule
        invalid_mol = Molecule(_rdkit_mol=nothing, valid=false, source="invalid")
        result = screen_molecule(invalid_mol, collection, molecule_id="invalid")

        @test result.molecule_id == "invalid"
        @test result.total_alerts == 0
        @test isempty(result.alert_matches)
        @test result.high_priority_count == 0
        @test result.passed == false  # Invalid molecules don't pass
    end

    @testset "screen_molecules - Vector Input" begin
        collection, molecules = setup_test_data()
        mol_ids = ["ethanol", "acetone", "benzene", "ethylamine", "propane", "acetophenone"]

        results = screen_molecules(molecules, collection, molecule_ids=mol_ids)

        @test length(results) == 6
        @test results[1].molecule_id == "ethanol"
        @test results[5].molecule_id == "propane"

        # Check specific results
        @test results[1].total_alerts == 1  # Ethanol
        @test results[2].total_alerts == 1  # Acetone
        @test results[3].total_alerts == 1  # Benzene
        @test results[4].total_alerts == 1  # Ethylamine
        @test results[5].total_alerts == 0  # Propane
        @test results[6].total_alerts == 2  # Acetophenone

        # Check pass/fail based on priority threshold (default 7)
        @test results[1].passed == false  # Alcohol priority 8
        @test results[2].passed == false  # Ketone priority 7
        @test results[3].passed == true   # Benzene priority 6
        @test results[4].passed == true   # Amine priority 5
        @test results[5].passed == true   # No alerts
        @test results[6].passed == false  # Ketone priority 7
    end

    @testset "screen_molecules - SMILES Vector Input" begin
        collection, _ = setup_test_data()
        smiles_list = ["CCO", "CC(=O)C", "CCC"]

        results = screen_molecules(smiles_list, collection)

        @test length(results) == 3
        @test results[1].total_alerts == 1  # Ethanol
        @test results[2].total_alerts == 1  # Acetone
        @test results[3].total_alerts == 0  # Propane
    end

    @testset "screen_molecules - Auto-generated IDs" begin
        collection, molecules = setup_test_data()

        results = screen_molecules(molecules[1:3], collection)  # No molecule_ids provided

        @test length(results) == 3
        @test results[1].molecule_id == "mol_1"
        @test results[2].molecule_id == "mol_2"
        @test results[3].molecule_id == "mol_3"
    end

    @testset "screen_molecules - ID Length Mismatch" begin
        collection, molecules = setup_test_data()

        @test_throws ArgumentError screen_molecules(
            molecules[1:3], collection, molecule_ids=["mol1", "mol2"]  # Wrong length
        )
    end

    @testset "get_alert_matches" begin
        collection, molecules = setup_test_data()

        result = screen_molecule(molecules[6], collection, molecule_id="acetophenone")  # Has ketone and benzene
        matched_alerts = get_alert_matches(result, collection)

        @test length(matched_alerts) == 2
        descriptions = [alert.description for alert in matched_alerts]
        @test "Ketone" in descriptions
        @test "Benzene" in descriptions
    end

    @testset "has_any_alert" begin
        collection, molecules = setup_test_data()

        result_with_alerts = screen_molecule(molecules[1], collection)  # Ethanol
        result_no_alerts = screen_molecule(molecules[5], collection)    # Propane

        @test has_any_alert(result_with_alerts) == true
        @test has_any_alert(result_no_alerts) == false
    end

    @testset "get_screening_summary" begin
        collection, molecules = setup_test_data()

        results = screen_molecules(molecules, collection, priority_threshold=7)
        summary = get_screening_summary(results)

        @test summary["total_molecules"] == 6
        @test summary["passed"] == 3  # Benzene, Amine, Propane pass (priority < 7 or no alerts)
        @test summary["failed"] == 3  # Ethanol, Acetone, Acetophenone fail (have priority >= 7)
        @test summary["pass_rate"] == 0.5
        @test summary["total_alert_instances"] == 6  # Total alerts across all molecules
        @test summary["max_alerts"] == 2  # Acetophenone has 2 alerts
        @test summary["avg_alerts_per_molecule"] == 1.0  # 6 alerts / 6 molecules
    end

    @testset "get_screening_summary - Empty Results" begin
        summary = get_screening_summary(ScreeningResult[])

        @test summary["total_molecules"] == 0
        @test summary["passed"] == 0
        @test summary["failed"] == 0
        @test summary["pass_rate"] == 0.0
        @test summary["avg_alerts_per_molecule"] == 0.0
        @test summary["max_alerts"] == 0
        @test summary["total_alert_instances"] == 0
    end

    @testset "Priority Threshold Effects" begin
        collection, molecules = setup_test_data()

        # Test with different priority thresholds
        results_threshold_8 = screen_molecules(molecules[1:2], collection, priority_threshold=8)
        results_threshold_6 = screen_molecules(molecules[1:2], collection, priority_threshold=6)

        # Ethanol has alcohol (priority 8), Acetone has ketone (priority 7)
        @test results_threshold_8[1].passed == false  # Alcohol >= 8
        @test results_threshold_8[2].passed == true   # Ketone < 8

        @test results_threshold_6[1].passed == false  # Alcohol >= 6
        @test results_threshold_6[2].passed == false  # Ketone >= 6
    end
end
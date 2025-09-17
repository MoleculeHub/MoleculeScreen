@testset "Alert Loading Tests" begin

    @testset "StructuralAlert Construction" begin
        alert = StructuralAlert(1, 1, "Test alert", "CCO", "TestSet", 8, 0)
        @test alert.rule_id == 1
        @test alert.rule_set == 1
        @test alert.description == "Test alert"
        @test alert.smarts == "CCO"
        @test alert.rule_set_name == "TestSet"
        @test alert.priority == 8
        @test alert.max == 0
    end

    @testset "AlertCollection Construction" begin
        alerts = [
            StructuralAlert(1, 1, "Alcohol", "[OH]", "Basic", 5, 0),
            StructuralAlert(2, 1, "Ketone", "C=O", "Basic", 6, 0),
            StructuralAlert(3, 2, "Benzene", "c1ccccc1", "Aromatic", 7, 0)
        ]

        collection = AlertCollection(alerts)
        @test length(collection.alerts) == 3
        @test length(collection.alert_dict) == 3
        @test haskey(collection.alert_dict, 1)
        @test haskey(collection.alert_dict, 2)
        @test haskey(collection.alert_dict, 3)
        @test collection.alert_dict[1].description == "Alcohol"
    end

    @testset "Load Alert Collection from DataFrame" begin
        # Create test DataFrame
        df = DataFrame(
            rule_id = [1, 2, 3],
            rule_set = [1, 1, 2],
            description = ["Test Alcohol", "Test Ketone", "Test Benzene"],
            smarts = ["[OH]", "C=O", "c1ccccc1"],
            rule_set_name = ["Basic", "Basic", "Aromatic"],
            priority = [5, 6, 7],
            max = [0, 0, 0]
        )

        collection = load_alert_collection(df)
        @test length(collection.alerts) == 3
        @test collection.alerts[1].description == "Test Alcohol"
        @test collection.alerts[2].smarts == "C=O"
        @test collection.alerts[3].rule_set_name == "Aromatic"
    end

    @testset "Invalid SMARTS Handling" begin
        # Create DataFrame with invalid SMARTS
        df = DataFrame(
            rule_id = [1, 2, 3],
            rule_set = [1, 1, 1],
            description = ["Valid", "Invalid", "Valid2"],
            smarts = ["[OH]", "invalid_smarts_pattern", "C=O"],
            rule_set_name = ["Basic", "Basic", "Basic"],
            priority = [5, 6, 7],
            max = [0, 0, 0]
        )

        # This should load successfully but skip the invalid SMARTS
        collection = load_alert_collection(df)
        @test length(collection.alerts) == 2  # Only valid patterns loaded
        @test collection.alerts[1].description == "Valid"
        @test collection.alerts[2].description == "Valid2"
    end

    @testset "get_alert_by_id" begin
        alerts = [
            StructuralAlert(1, 1, "Alcohol", "[OH]", "Basic", 5, 0),
            StructuralAlert(2, 1, "Ketone", "C=O", "Basic", 6, 0)
        ]
        collection = AlertCollection(alerts)

        alert1 = get_alert_by_id(collection, 1)
        @test alert1 !== nothing
        @test alert1.description == "Alcohol"

        alert_missing = get_alert_by_id(collection, 999)
        @test alert_missing === nothing
    end

    @testset "filter_alerts_by_rule_set" begin
        alerts = [
            StructuralAlert(1, 1, "Alcohol", "[OH]", "Basic", 5, 0),
            StructuralAlert(2, 1, "Ketone", "C=O", "Basic", 6, 0),
            StructuralAlert(3, 2, "Benzene", "c1ccccc1", "Aromatic", 7, 0)
        ]
        collection = AlertCollection(alerts)

        basic_alerts = filter_alerts_by_rule_set(collection, "Basic")
        @test length(basic_alerts) == 2
        @test all(alert.rule_set_name == "Basic" for alert in basic_alerts)

        aromatic_alerts = filter_alerts_by_rule_set(collection, "Aromatic")
        @test length(aromatic_alerts) == 1
        @test aromatic_alerts[1].description == "Benzene"

        missing_alerts = filter_alerts_by_rule_set(collection, "NonExistent")
        @test length(missing_alerts) == 0
    end

    @testset "filter_alerts_by_priority" begin
        alerts = [
            StructuralAlert(1, 1, "Low", "[OH]", "Basic", 3, 0),
            StructuralAlert(2, 1, "Medium", "C=O", "Basic", 6, 0),
            StructuralAlert(3, 1, "High", "c1ccccc1", "Basic", 8, 0)
        ]
        collection = AlertCollection(alerts)

        high_priority = filter_alerts_by_priority(collection, 7)
        @test length(high_priority) == 1
        @test high_priority[1].description == "High"

        medium_plus = filter_alerts_by_priority(collection, 6)
        @test length(medium_plus) == 2
        @test all(alert.priority >= 6 for alert in medium_plus)

        all_alerts = filter_alerts_by_priority(collection, 1)
        @test length(all_alerts) == 3
    end

    @testset "get_rule_set_names" begin
        alerts = [
            StructuralAlert(1, 1, "Test1", "[OH]", "Basic", 5, 0),
            StructuralAlert(2, 1, "Test2", "C=O", "Basic", 6, 0),
            StructuralAlert(3, 2, "Test3", "c1ccccc1", "Aromatic", 7, 0),
            StructuralAlert(4, 3, "Test4", "[N]", "Nitrogen", 8, 0)
        ]
        collection = AlertCollection(alerts)

        rule_sets = get_rule_set_names(collection)
        @test length(rule_sets) == 3
        @test "Basic" in rule_sets
        @test "Aromatic" in rule_sets
        @test "Nitrogen" in rule_sets
        @test issorted(rule_sets)  # Should be sorted
    end

    @testset "get_priority_distribution" begin
        alerts = [
            StructuralAlert(1, 1, "Test1", "[OH]", "Basic", 5, 0),
            StructuralAlert(2, 1, "Test2", "C=O", "Basic", 5, 0),
            StructuralAlert(3, 1, "Test3", "c1ccccc1", "Basic", 7, 0),
            StructuralAlert(4, 1, "Test4", "[N]", "Basic", 8, 0),
            StructuralAlert(5, 1, "Test5", "[S]", "Basic", 8, 0)
        ]
        collection = AlertCollection(alerts)

        dist = get_priority_distribution(collection)
        @test dist[5] == 2
        @test dist[7] == 1
        @test dist[8] == 2
        @test length(dist) == 3
    end

    @testset "Missing Columns Error" begin
        # DataFrame missing required columns
        incomplete_df = DataFrame(
            rule_id = [1, 2],
            description = ["Test1", "Test2"]
            # Missing other required columns
        )

        @test_throws ArgumentError load_alert_collection(incomplete_df)
    end

    @testset "File Not Found Error" begin
        @test_throws ArgumentError load_alert_collection("nonexistent_file.csv")
    end
end
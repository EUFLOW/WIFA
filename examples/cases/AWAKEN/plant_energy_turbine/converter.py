import os
import re
import sys

import yaml


class NoCurlyBracesDumper(yaml.SafeDumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(NoCurlyBracesDumper, self).increase_indent(
            flow=False, indentless=False
        )


def extract_array_from_text(text, array_name):
    # Regex to find the np.array assignment and extract the array content
    pattern = rf"{array_name}\s*=\s*np\.array\(\[\[(.*?)\]\]\)"
    match = re.search(pattern, text, re.DOTALL)

    if not match:
        raise ValueError(f"Could not find array {array_name} in the file")

    # Extract the array string and clean it up
    array_data = match.group(1).replace("\n", "").replace(" ", "")

    # Split and convert to float
    array_data = [
        [float(num) for num in row.split(",")] for row in array_data.split("],[")
    ]

    return array_data


def find_variable_prefix(text):
    # Regex to dynamically find any power curve and ct curve variable names
    power_curve_match = re.search(r"(\w+)_power_curve\s*=\s*np\.array", text)
    ct_curve_match = re.search(r"(\w+)_ct_curve\s*=\s*np\.array", text)

    if power_curve_match and ct_curve_match:
        prefix = power_curve_match.group(1)
        return prefix
    else:
        raise ValueError("Could not find power_curve or ct_curve variables in the file")


def process_turbine_file(file_path):
    # Read the entire Python file as a text file
    with open(file_path, "r") as f:
        file_content = f.read()

    # Find the variable prefix dynamically
    prefix = find_variable_prefix(file_content)

    # Extract the power and Ct curves
    try:
        power_curve = extract_array_from_text(file_content, f"{prefix}_power_curve")
        ct_curve = extract_array_from_text(file_content, f"{prefix}_ct_curve")
    except ValueError as e:
        print(f"Error: {e}")
        return

    # Extract power values and wind speeds from power curve
    power_wind_speeds = [row[0] for row in power_curve]
    power_values = [row[1] for row in power_curve]

    # Extract Ct values and wind speeds from Ct curve
    ct_wind_speeds = [row[0] for row in ct_curve]
    ct_values = [row[1] for row in ct_curve]

    # Define the data structure for YAML with list formatting in brackets
    turbine_data = {
        "name": prefix.upper().replace("P", ".").replace("MW", " MW"),
        "performance": {
            "power_values": power_values,
            "power_wind_speeds": power_wind_speeds,
            "Ct_curve": {"Ct_values": ct_values, "Ct_wind_speeds": ct_wind_speeds},
        },
        "cutin_wind_speed": 3.0,
        "cutout_wind_speed": 25.0,
        "hub_height": 80.0,  # Modify this value based on turbine data, if needed
        "rotor_diameter": 103.0,  # Modify this value based on turbine data, if needed
    }

    # Convert to YAML format with lists in brackets but no curly braces for dictionaries
    yaml_string = yaml.dump(
        turbine_data,
        Dumper=NoCurlyBracesDumper,
        sort_keys=False,
        default_flow_style=None,
    )
    print(f"\nYAML Output for {prefix}:\n")
    print(yaml_string)

    # Save to a YAML file
    output_file = f"{prefix}.yaml".strip("mw")
    with open(output_file, "w") as file:
        yaml.dump(
            turbine_data,
            file,
            Dumper=NoCurlyBracesDumper,
            sort_keys=False,
            default_flow_style=None,
        )
    print(f"Saved {prefix} data to {output_file}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python converter.py <list_of_turbine_files>")
        sys.exit(1)

    # Iterate over each turbine file provided as input
    for turbine_file in sys.argv[1:]:
        print(f"Processing turbine file: {turbine_file}")
        process_turbine_file(turbine_file)


if __name__ == "__main__":
    main()

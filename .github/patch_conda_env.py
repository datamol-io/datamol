import sys
import argparse
import yaml


def main(env_path, new_deps):

    # Process dependencies
    new_deps = [dep.split("=") for dep in new_deps]

    # Load the conda env file
    with open(env_path) as f:
        env_spec = yaml.load(f, Loader=yaml.SafeLoader)

    deps = env_spec["dependencies"]
    for name, version in new_deps:

        # Find whether the package already exists
        existing = list(
            filter(
                lambda x: True if isinstance(x, str) and x.startswith(name) else False,
                deps,
            )
        )

        # Remove the existing package(s)
        [deps.pop(deps.index(e)) for e in existing]

        # Add the new package spec if the spec is not None
        if version != "DELETE":
            if version == "":
                deps.append(name)
            else:
                deps.append(f"{name} ={version}")

    # Add the new dep list to the eocnda env spec
    env_spec["dependencies"] = deps

    # Serialize back to YAML and print on stdout
    sys.stdout.write(yaml.dump(env_spec))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Patch a conda env file")
    parser.add_argument(
        "--env",
        metavar="CONDA_ENV_PATH",
        type=str,
        help="Path to conda env file.",
    )
    parser.add_argument(
        "-d",
        "--deps",
        nargs="+",
        help="Dependencies to patch.",
    )

    args = parser.parse_args()
    main(args.env, args.deps)

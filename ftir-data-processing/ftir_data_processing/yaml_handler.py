import yaml


def write_yaml(filepath, metadata):
    with open(filepath, "w", encoding="utf-8") as f:
        yaml.safe_dump(metadata, f, sort_keys=False, allow_unicode=True)


def read_yaml(filepath: str) -> yaml:
    with open(filepath, 'r') as f:
        metadata = yaml.safe_load(f)
    return metadata


if __name__ == "__main__":
    yaml_fn = 'D:\\Promotie\\scripts\\ftir-data-processing\\templates\\IRRAS\\NH3_plasma.yaml'
    reference_metadata = read_yaml(yaml_fn)

    print(reference_metadata)

import requests
import base64
import json
import sys

from ord_schema.proto import reaction_pb2
from rdkit import Chem
from google.protobuf.json_format import MessageToDict


# =========================
# SMARTS PATTERNS
# =========================
AMINE_SMARTS = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
ARYL_HALIDE_SMARTS = Chem.MolFromSmarts("[c]~[F,Cl,Br,I]")
CARBOXYLIC_ACID_SMARTS = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

METAL_ATOMIC_NUMBERS = {
    3, 11, 12, 19, 20, 26, 27, 28, 29, 30,
    37, 38, 44, 45, 46, 47, 48,
    55, 56, 76, 77, 78, 79, 80
}


# =========================
# HELPER FUNCTIONS
# =========================
def get_smiles(component):
    for identifier in component.identifiers:
        if identifier.type == reaction_pb2.CompoundIdentifier.SMILES:
            return identifier.value
    return None


def safe_mol_from_smiles(smiles):
    try:
        return Chem.MolFromSmiles(smiles) if smiles else None
    except Exception:
        return None


def is_metal(mol):
    if mol is None:
        return False
    return any(atom.GetAtomicNum() in METAL_ATOMIC_NUMBERS for atom in mol.GetAtoms())


def get_best_identifier(component):
    for ident in component.identifiers:
        if ident.type == reaction_pb2.CompoundIdentifier.SMILES:
            return "SMILES", ident.value
    for ident in component.identifiers:
        if ident.type == reaction_pb2.CompoundIdentifier.NAME:
            return "NAME", ident.value
    return "UNKNOWN", "Unknown"


# =========================
# CLASSIFICATION LOGIC
# =========================
def classify_component(component, role, input_key=""):
    smiles = get_smiles(component)
    mol = safe_mol_from_smiles(smiles)
    key = input_key.lower()

    categories = set()

    # 1. Role override
    if role == reaction_pb2.ReactionRole.SOLVENT:
        return {"Solvent"}

    # 2. Input key hints
    if "carboxylic acid" in key or "acid" in key:
        categories.add("Carboxylic Acid")

    if "amine" in key:
        categories.add("Amine")

    if "base" in key:
        categories.add("Base")

    if "additive" in key:
        categories.add("Additive")

    if "activation" in key or "coupling agent" in key:
        categories.add("Activation Agent")

    if "ligand" in key:
        categories.add("Ligand")

    if "catalyst" in key or "metal" in key:
        categories.add("Metal" if is_metal(mol) else "Ligand")

    # 3. Role-based inference
    if role == reaction_pb2.ReactionRole.CATALYST and not categories:
        categories.add("Metal" if is_metal(mol) else "Ligand")

    if role == reaction_pb2.ReactionRole.REACTANT and mol:
        if mol.HasSubstructMatch(AMINE_SMARTS):
            categories.add("Amine")
        if mol.HasSubstructMatch(ARYL_HALIDE_SMARTS):
            categories.add("Aryl Halide")
        if mol.HasSubstructMatch(CARBOXYLIC_ACID_SMARTS):
            categories.add("Carboxylic Acid")

    # 4. Generic reagent fallback
    if role == reaction_pb2.ReactionRole.REAGENT and not categories:
        categories.add("Base")

    # 5. M-series handling
    if key.startswith("m"):
        for part in key.split("_"):
            if part.startswith("m") and part[1:].isdigit():
                categories.add(part.upper())

    return categories


# =========================
# MAIN SCRIPT
# =========================
def main():
    datasets_to_process = []

    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            if "ord_dataset-" in arg:
                idx = arg.find("ord_dataset-")
                dataset_id = arg[idx:].split("/")[0]
                datasets_to_process.append({"dataset_id": dataset_id})

    if not datasets_to_process:
        datasets_url = "https://open-reaction-database.org/api/datasets"
        try:
            datasets_to_process = requests.get(datasets_url).json()[:2]
        except Exception as e:
            print(f"Failed to fetch datasets: {e}")
            return

    organized_data = {
        "Base": [],
        "Solvent": [],
        "Amine": [],
        "Aryl Halide": [],
        "Metal": [],
        "Ligand": [],
        "Carboxylic Acid": [],
        "Additive": [],
        "Activation Agent": [],
        "M1": [], "M2": [], "M3": [], "M4": [],
        "M5": [], "M6": [], "M7": [], "M8": [], "M9": []
    }

    for dataset in datasets_to_process:
        dataset_id = dataset["dataset_id"]
        print(f"Processing dataset {dataset_id}")

        query_url = "https://open-reaction-database.org/api/query"
        params = {"dataset_id": dataset_id, "limit": 50}

        try:
            results = requests.get(query_url, params=params).json()
        except Exception as e:
            print(f"Query failed for {dataset_id}: {e}")
            continue

        for result in results:
            reaction_id = result["reaction_id"]

            try:
                proto_bytes = base64.b64decode(result["proto"])
                reaction = reaction_pb2.Reaction.FromString(proto_bytes)
            except Exception as e:
                print(f"Failed to decode reaction {reaction_id}: {e}")
                continue

            for input_key, input_val in reaction.inputs.items():
                for component in input_val.components:
                    raw = MessageToDict(
                        component,
                        preserving_proto_field_name=True,
                        use_integers_for_enums=False
                    )

                    role_name = raw.get("reaction_role", "UNSPECIFIED")
                    categories = classify_component(
                        component,
                        component.reaction_role,
                        input_key
                    )

                    identifier_type, value = get_best_identifier(component)

                    for cat in categories:
                        organized_data.setdefault(cat, []).append({
                            "reaction_id": reaction_id,
                            "input_key": input_key,
                            "reaction_role": role_name,
                            "identifier_type": identifier_type,
                            "value": value
                        })

    output = {
        "metadata": {
            "source": "Open Reaction Database",
            "total_components": sum(len(v) for v in organized_data.values())
        },
        "classified_components": organized_data
    }

    with open("ord_data.json", "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)

    print("Saved ord_data.json")


if __name__ == "__main__":
    main()

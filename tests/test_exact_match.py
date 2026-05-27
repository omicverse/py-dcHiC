"""The parity gate, as a pytest test.

Runs the dcHiC R reference and the pydchic candidate on the same fixture, applies the
pre-registered class-aware metrics from the manifest, asserts both gates pass.

Do NOT loosen a threshold to make this pass — fix the candidate.
"""

import json
import os
import subprocess
from pathlib import Path

import pytest
import yaml

from .parity_metrics import compute_parity, is_pass

HERE = Path(__file__).parent
PORT_DIR = HERE.parent
MANIFEST_PATH = PORT_DIR / "data" / "manifest.yaml"
R_CONDA_ENV = os.environ.get("R_TEST_ENV", "/home/shengmao/.local/share/mamba/envs/dchic")
PY_CONDA_ENV = os.environ.get("PYTHON_TEST_ENV", "/home/shengmao/miniconda/envs/rebuild-py")


@pytest.fixture(scope="session")
def manifest():
    with open(MANIFEST_PATH) as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="session")
def fixture_path(manifest):
    path = PORT_DIR / manifest["fixture"]["path"]
    if not path.exists():
        subprocess.run(
            [f"{PY_CONDA_ENV}/bin/python", str(HERE / "build_fixture.py")], check=True
        )
    return path


@pytest.fixture(scope="session")
def reference_output(manifest, fixture_path):
    cache = PORT_DIR / "data" / "reference_output.json"
    if cache.exists():
        return json.loads(cache.read_text())
    ref_script = PORT_DIR / manifest["reference_command"]
    subprocess.run(
        [f"{R_CONDA_ENV}/bin/Rscript", str(ref_script), str(fixture_path), str(cache)],
        check=True,
    )
    return json.loads(cache.read_text())


@pytest.fixture(scope="session")
def candidate_output(manifest, fixture_path):
    cache = PORT_DIR / "data" / "candidate_output.json"
    runner = HERE / "_run_candidate.py"
    env = dict(os.environ, PYTHONPATH=str(PORT_DIR))
    subprocess.run(
        [f"{PY_CONDA_ENV}/bin/python", str(runner), str(fixture_path), str(cache)],
        check=True, env=env,
    )
    return json.loads(cache.read_text())


def test_parity_against_R(manifest, reference_output, candidate_output):
    failures, report = [], []
    for spec in manifest["outputs"]:
        name = spec["name"]
        ref_key = spec["location_reference"].lstrip("$.")
        cand_key = spec["location_candidate"]
        cls = spec.get("metric", manifest["algorithm_class"])
        threshold = spec.get("threshold", manifest["parity_threshold"])
        metric = compute_parity(reference_output[ref_key], candidate_output[cand_key], cls)
        passed = is_pass(metric, cls, threshold)
        report.append(f"[parity] {name} ({cls}): {metric} threshold={threshold} pass={passed}")
        if not passed:
            failures.append((name, metric, threshold))
    print("\n".join(report))
    assert not failures, f"Parity gate failed: {failures}"

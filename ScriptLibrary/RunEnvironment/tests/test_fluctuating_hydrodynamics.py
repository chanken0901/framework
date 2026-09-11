import copy
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from case_input import CaseInputError, render_nse, render_nse_multicomponent
from case_configuration import resolve_case_configuration
import test_case_input


class FluctuatingTests(unittest.TestCase):
    def setUp(self):
        from yaml_support import load_yaml
        root = Path(__file__).resolve().parents[3]
        self.manifest = load_yaml(root / "SolverLibrary/NSE/solver_manifest.yaml")
        self.case = test_case_input.NseCaseInputTests.case()
        self.case["physics"]["nse"]["reynolds_number"] = 100.0
        self.case["numerics"]["viscous_scheme"] = "central6"
        self.case["numerics"]["convective_scheme"] = "keep6"
        self.case["time"]["use_fixed_dt"] = True
        self.case["time"]["dt"] = 1.e-4
        self.case["physics"]["fluctuating_hydrodynamics"] = {
            "model": "landau_lifshitz", "boltzmann_number": 1.e-8, "seed": 42}

    def test_render(self):
        for profile in ("cpu_mpi", "cuda_single", "cuda_mpi"):
            text = render_nse(self.case, self.manifest, profile)
            self.assertIn("fh_enabled = .true.", text)
            self.assertIn("fh_seed = 42", text)

    def test_disabled_and_legacy(self):
        del self.case["physics"]["fluctuating_hydrodynamics"]
        self.assertNotIn("fh_enabled", render_nse(self.case, self.manifest, "cpu_mpi"))
        self.case["physics"]["fluctuating_hydrodynamics"] = {"enabled": False}
        self.assertIn("fh_enabled = .false.", render_nse(self.case, self.manifest, "cpu_mpi"))

    def test_standard_controls_are_preserved(self):
        from profile_selection import select_case_profile
        case = copy.deepcopy(self.case)
        case["forcing"] = {"type": "petersen_livescu", "petersen_livescu": {
            "spectrum": "low_wavenumber", "target_dissipation": 0.1,
            "dilatational_ratio": 0.0, "k_cutoff": 2.5}}
        case["flow"]["type"] = "hit"
        case.get("solver", {}).pop("profile", None)
        profile, _ = select_case_profile(case, self.manifest, {"model":"nse",
            "profile":"cpu_mpi", "available_profiles":["cpu_mpi","cpu_mpi_2decomp_fftw"]})
        self.assertEqual(profile, "cpu_mpi_2decomp_fftw")
        enabled = render_nse(case, self.manifest, profile)
        case["physics"]["fluctuating_hydrodynamics"] = {"enabled":False}
        disabled = render_nse(case, self.manifest, profile)
        del case["physics"]["fluctuating_hydrodynamics"]
        legacy = render_nse(case, self.manifest, profile)
        def without_fh(text):
            return '\n'.join(line for line in text.splitlines() if not line.strip().startswith('fh_'))
        self.assertEqual(without_fh(enabled), legacy.rstrip())
        self.assertEqual(without_fh(disabled), legacy.rstrip())

    def test_reject_ignored_controls(self):
        for key in ('chemistry','geometry','thermodynamics','transport'):
            case=copy.deepcopy(self.case); case[key]={}
            with self.subTest(key=key), self.assertRaisesRegex(CaseInputError,key):
                render_nse(case,self.manifest,'cpu_mpi')
        baseline=test_case_input.MulticomponentFoundationInputTests.reactive_case()
        for section, key, value in [('forcing','type','petersen_livescu'),
                ('time','t_max',1.0),('time','use_fixed_dt',True),
                ('time','output_frequency',5),('output','format','slf'),
                ('output','write_initial',True),('output','directory','output')]:
            case=copy.deepcopy(baseline); case.setdefault(section,{})[key]=value
            with self.subTest(key=key), self.assertRaises(CaseInputError):
                render_nse_multicomponent(case,'cpu_serial_reactive')

    def test_template_overlay_contract(self):
        from case_template_overlay import apply_template_overrides
        text='flow:\n  type: taylor_green # retained\nforcing:\n  type: none\n'
        out=apply_template_overrides(text,{'flow.type':'hit','schema_version':2})
        self.assertIn('# retained',out)
        self.assertIn('forcing:\n  type: none',out)
        for changes in ({'flow.typo':1},{'flow':{}},{'bad/path':1}):
            with self.assertRaises(ValueError): apply_template_overrides(text,changes)

    def test_invalid(self):
        for key, value in [("boltzmann_number", -1), ("boltzmann_number", float("nan")),
                           ("boltzmann_number", True), ("seed", -1), ("seed", 2**31),
                           ("seed", True), ("enabled", "yes"), ("model", "random_force"),
                           ("typo", 1)]:
            case = copy.deepcopy(self.case)
            case["physics"]["fluctuating_hydrodynamics"][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(CaseInputError):
                render_nse(case, self.manifest, "cpu_mpi")

    def test_unsupported(self):
        with self.assertRaises(CaseInputError):
            render_nse_multicomponent(self.case)
        for section, key, value in [("numerics", "viscous_scheme", "none")]:
            case = copy.deepcopy(self.case)
            case[section][key] = value
            with self.subTest(key=key), self.assertRaises(CaseInputError):
                render_nse(case, self.manifest, "cpu_mpi")

    def test_all_convection_and_time_modes(self):
        for scheme in ("keep2", "keep6", "weno5z_roe", "hybrid"):
            for fixed in (False, True):
                case = copy.deepcopy(self.case)
                case["numerics"]["convective_scheme"] = scheme
                case["time"]["use_fixed_dt"] = fixed
                with self.subTest(scheme=scheme, fixed=fixed):
                    self.assertIn("fh_enabled = .true.", render_nse(case, self.manifest, "cpu_mpi"))

    def test_sidecar(self):
        import json
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case = copy.deepcopy(self.case)
            config = case["physics"].pop("fluctuating_hydrodynamics")
            case.update(schema_version=2, extensions={"fluctuating_hydrodynamics": "fh.yaml"})
            (root / "case.yaml").write_text(json.dumps(case))
            (root / "fh.yaml").write_text(json.dumps({"schema_version": 1,
                "extension": "fluctuating_hydrodynamics", "config": config}))
            result = resolve_case_configuration(root / "case.yaml")
            self.assertEqual(result.document["physics"]["fluctuating_hydrodynamics"], config)
            self.assertIn("fh_seed = 42", render_nse(result.document, self.manifest, "cpu_mpi"))

    def test_generated_cuda_environments(self):
        from argparse import Namespace
        from unittest.mock import patch
        from prepare_environment import prepare
        from yaml_support import load_yaml
        script = Path(__file__).resolve().parents[1]
        design_path = script / "environment.nse_fluctuating.cuda.yaml"
        for use_mpi in (False, True):
            with self.subTest(use_mpi=use_mpi), tempfile.TemporaryDirectory() as directory:
                design = load_yaml(design_path)
                design["parallel"]["use_mpi"] = use_mpi
                def load_with_design(path):
                    return copy.deepcopy(design) if Path(path).resolve() == design_path.resolve() else load_yaml(path)
                args = Namespace(design=str(design_path), framework_root=str(script.parents[1]),
                    output=str(Path(directory) / "llns"), model=None, profile=None,
                    case_id=None, overwrite=False, archive=False, archive_format=None, dry_run=False)
                with patch("prepare_environment.sync_environment_case"), \
                     patch("prepare_environment.load_yaml", side_effect=load_with_design):
                    generated = prepare(args)
                source = generated / "SolverLibrary/NSE"
                self.assertTrue((source / "src/gpu/nse_cuda_fluctuating.cuh").is_file())
                self.assertTrue((source / "src/extensions/fluctuating/mod_nse_fluctuating.f90").is_file())
                self.assertEqual((source / "tests/test_cuda_fluctuating_mpi.f90").is_file(), use_mpi)
                config = resolve_case_configuration(generated / "cases/case0001/case.yaml")
                self.assertIn("fluctuating_hydrodynamics", config.extension_paths)
                text = (generated / "cases/case0001/input.dat").read_text()
                self.assertIn("fh_enabled = .true.", text)

    def test_generated_environment(self):
        from argparse import Namespace
        from unittest.mock import patch
        from prepare_environment import prepare
        script = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            args = Namespace(design=str(script / "environment.nse_fluctuating.yaml"),
                framework_root=str(script.parents[1]), output=str(Path(directory) / "llns"),
                model=None, profile=None, case_id=None, overwrite=False, archive=False,
                archive_format=None, dry_run=False)
            # This integration test must not touch the user's global case index.
            with patch("prepare_environment.sync_environment_case"):
                generated = prepare(args)
            source = generated / "SolverLibrary/NSE"
            self.assertTrue((source / "src/extensions/fluctuating/mod_fh_random.f90").is_file())
            self.assertTrue((source / "src/extensions/fluctuating/mod_nse_fluctuating.f90").is_file())
            config = resolve_case_configuration(generated / "cases/case0001/case.yaml")
            self.assertIn("fluctuating_hydrodynamics", config.extension_paths)
            self.assertIn("hit", config.document["flow"])
            self.assertIn("imported_turbulence", config.document["flow"])
            self.assertEqual(config.document["forcing"]["type"], "none")
            self.assertIn("hybrid", config.document["numerics"])
            template_text=(generated / "templates/case_template.yaml").read_text(encoding="utf-8")
            self.assertIn("Petersen-Livescu",template_text)
            self.assertIn("fh_enabled = .true.",
                (generated / "cases/case0001/input.dat").read_text())
            # Optional toolchain-dependent integration check. All artifacts stay
            # inside this temporary directory; no user ResearchRuns writes.
            import os
            if os.environ.get("NSE_FH_TEST_BUILD") == "1":
                import subprocess
                binary = generated / "build-check"
                command = ["cmake", "-S", str(source), "-B", str(binary), "-G", "Ninja",
                    "-DNSE_USE_MPI=ON", "-DNSE_ENABLE_OPENMP=ON", "-DNSE_GPU_BACKEND=none",
                    "-DNSE_INIT_FFT_BACKEND=none", "-DNSE_FORCING_FFT_BACKEND=none",
                    "-DNSE_VISCOUS_SCHEME=central6", "-DBUILD_TESTING=ON"]
                if os.name == "nt":
                    command.append("-DNSE_MPI_PROVIDER=MSMPI")
                subprocess.run(command, check=True)
                subprocess.run(["cmake", "--build", str(binary), "-j", "4"], check=True)
                subprocess.run(["ctest", "--test-dir", str(binary), "-R", "fluctuating|fh_reject",
                    "--output-on-failure"], check=True)

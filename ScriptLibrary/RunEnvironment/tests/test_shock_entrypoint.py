from pathlib import Path
import sys
import tempfile
import shutil
import json
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[3]
SCRIPT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SCRIPT))
sys.path.insert(0,str(ROOT/'SolverLibrary/NSE/tests'))
import analyze_shock_turbulence as entry
import postprocess_case
from test_shock_analysis import write_slf,state


class EntryTests(unittest.TestCase):
    def test_dedicated_entrypoint_and_manifest(self):
        root_source=SCRIPT.parents[1]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);output=root/'cases/case1/output';output.mkdir(parents=True)
            tools=root/'SolverLibrary/NSE/tools';tools.mkdir(parents=True)
            for name in ('nse_shock_analysis.py','nse_turbulence_statistics.py','slf_to_paraview_merged_cropghost.py'):
                shutil.copy2(root_source/'SolverLibrary/NSE/tools'/name,tools/name)
            write_slf(output/'field_000010.slf',state())
            meta=output/'meta.json';meta.write_text(json.dumps(dict(grid=[40,4,4],spacing=[1,1,1])))
            context=dict(root=root,model='nse',case_root=output.parent,input_dir=output,meta=meta,gamma=1.4)
            (root/'environment.lock.json').write_text(json.dumps(dict(model='nse',case_directory='cases/case1')))
            (output.parent/'case.yaml').write_text('physics:\n  nse:\n    gamma: 1.4\n')
            args=['--direction','1','--search','10','30','--upstream','2','8','--downstream','2','8']
            with patch.object(postprocess_case,'__file__',str(root/'tools/postprocess_case.py')):
                self.assertEqual(entry.main(args+['--dry-run']),0)
                self.assertFalse((output.parent/'shock_analysis').exists())
                self.assertEqual(entry.main(args),0)
                self.assertTrue((output.parent/'shock_analysis/analysis.json').is_file())
            with patch.object(entry,'_case_context',return_value=context):
                context['model']='gpe'
                self.assertEqual(entry.main(args),1)
        self.assertIn('"analyze_shock_turbulence.py"',(SCRIPT/'prepare_environment.py').read_text(encoding='utf-8'))
        self.assertIn('tools/nse_shock_analysis.py',(root_source/'SolverLibrary/NSE/solver_manifest.yaml').read_text(encoding='utf-8'))

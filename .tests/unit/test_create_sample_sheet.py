import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_create_sample_sheet():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/create_sample_sheet/data")
        # expected_path = PurePosixPath(".tests/unit/python_script_test/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/create_sample_sheet.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)
        
        # run function
        sys.path.insert(0, workdir)
        from create_sample_sheet import main
        main(
            isolate_sheets=os.path.join(workdir, 'isolate_sheets.json'),
            fasta_prefix="",
            dirout=os.path.join(workdir),
        )

        # Insert your tests here
        assert os.path.exists(os.path.join(workdir, "Salmonella_enterica.tsv"))
        assert os.path.exists(os.path.join(workdir, "Campylobacter_spp.tsv"))

import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_make_isolate_sheet():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/make_isolate_sheet/data")
        expected_path = PurePosixPath(".tests/unit/make_isolate_sheet/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/make_isolate_sheet.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from make_isolate_sheet import main
        main(
            validation=os.path.join(workdir, 'validation.json'),
            calling=os.path.join(workdir, 'calling.json'),
            charak=os.path.join(workdir, 'charak.json'),
            pathout=os.path.join(workdir, "result.json"),
        )

        # Insert your tests here
        with open(
            os.path.join(workdir, 'result.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'merged.json'), 'r'
        ) as expect:
            assert load(res) == load(expect)

dir="$(dirname "$0")"
python "$dir/cleanup_db.py"
python "$dir/make_test_db.py"
echo "DONE"
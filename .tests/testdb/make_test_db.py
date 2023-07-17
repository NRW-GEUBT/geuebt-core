import os
import pathlib
import glob
import json
from pymongo import MongoClient


HOST = 'localhost'
PORT = 27017
ISOLATES = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data/isolate_sheets")

CLUSTERS = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data/clusters")


def listjson(folder):
    return list(pathlib.Path(folder).glob('*.json'))


def clean_load(filepath):
    with open(filepath, 'r') as fi:
        d = json.load(fi)
    return d


def paths_to_json(pathlist):
    return  [clean_load(filepath) for filepath in pathlist]


# Set up collections
client = MongoClient(HOST, PORT)
db = client.testing_db
isolates_coll = db.isolates
clusters_coll = db.clusters


# Insert documents
_ = isolates_coll.insert_many(paths_to_json(listjson(ISOLATES)))
_ = clusters_coll.insert_many(paths_to_json(listjson(CLUSTERS)))

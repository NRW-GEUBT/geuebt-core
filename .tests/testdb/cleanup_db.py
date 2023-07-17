from pymongo import MongoClient


HOST = 'localhost'
PORT = 27017


client = MongoClient(HOST, PORT)
db = client.testing_db
isolates_coll = db.isolates
clusters_coll = db.clusters


isolates_coll.drop()
clusters_coll.drop()

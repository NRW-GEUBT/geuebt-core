from pymongo import MongoClient


HOST = 'localhost'
PORT = 27017


def main(host, port):
    client = MongoClient(host, port)
    db = client.testing_db
    isolates_coll = db.isolates
    clusters_coll = db.clusters

    isolates_coll.drop()
    clusters_coll.drop()


if __name__ == '__main__':
    main(HOST, PORT)
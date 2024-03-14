from pymongo import MongoClient


HOST = 'localhost'
PORT = 27017
DATABASE = 'geuebt-test'


def main(host, port):
    client = MongoClient(host, port)
    db = client[DATABASE]
    db.isolates.drop()
    db.clusters.drop()
    db.runs.drop()
    db.sequences.drop()


if __name__ == '__main__':
    main(HOST, PORT)

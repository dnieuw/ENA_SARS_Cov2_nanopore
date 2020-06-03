#!/usr/bin/env python

import sys
import datetime

from pymongo import MongoClient


def main():
    sample = DB.samples.find_one({'id': RUN})
    sample[UPDATE_FIELD]['date'].append(
        datetime.datetime.now().strftime("%d %B, %Y - %H:%M:%S"))
    sample[UPDATE_FIELD]['status'].append(STATUS)
    DB.samples.update_one({'id': RUN}, {'$set': sample})


if __name__ == "__main__":
    # Get run name from command line
    RUN = sys.argv[1]
    if '_' in RUN:
        RUN = RUN.split('_')[0]

    # Get status from command line
    STATUS = sys.argv[2]
    UPDATE_FIELD = sys.argv[3]

    # Getting access to MongoDB
    CLIENT = MongoClient('mongodb://samples-logs-db-svc')
    DB = CLIENT.samples
    main()

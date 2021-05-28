import os
import sys

import pymongo
import pymongo.errors


def dbsnp_calculator(human_file, mongodb_info):
    """
    Annotate a .human file with dbSNP is, REF and ALT allele frequencies
    :param str human_file: The .human file to be annotated
    :param ("MongoDB_Info", "ip, port, dbsnp_ver") mongodb_info: A named tuple, containing (ip, port, dbsnp version)
    """

    client = pymongo.MongoClient(host=mongodb_info.ip, port=mongodb_info.port, connect=False)
    try:
        client.admin.command('ismaster')
        if mongodb_info.genome_version == 'hg19':
            dbsnp_database = client.dbsnp
            """:type: pymongo.database"""
        elif mongodb_info.genome_version == 'hg38':
            dbsnp_database = client.dbsnp_hg38
            """:type: pymongo.database"""
        dbsnp_collection_text = 'dbsnp_database.dbsnp{}'.format(mongodb_info.dbsnp_ver)
        dbsnp_collection = eval(dbsnp_collection_text)
        if dbsnp_collection.count() == 0:
            raise pymongo.errors.OperationFailure("Collection '{}' not present in the database 'dbsnp'".format(
                "dbsnp{}".format(mongodb_info.dbsnp_ver)))

        with open(human_file, "r") as hfile:
            header = hfile.readline().rstrip()
            newheader = header + "\tdbSNP ID\tFreq Ref (dbSNP%s)\tFreq Alt (dbSNP%s)\n" % \
                                 (mongodb_info.dbsnp_ver, mongodb_info.dbsnp_ver)

            outfile = os.path.splitext(os.path.abspath(human_file))[0] + "_freq.human"
            with open(outfile, "w") as out:
                print("--- Adding allelic frequencies from dbSNP ver. {} ---".format(mongodb_info.dbsnp_ver))
                out.write(newheader)

                for hline in hfile:
                    hline = hline.rstrip()  # rstrip removes pesky newline characters
                    htokens = hline.split("\t")
                    chr = htokens[0]
                    pos = int(htokens[1])
                    ref = htokens[3]
                    alt = htokens[4]
                    # print(htokens)


                    rsid = "."
                    rsid_list = []
                    """:type : list[str]"""
                    ref_freq = "."
                    alt_freq = "."
                    mongo_cursor = dbsnp_collection.find(
                        {"$and": [{"chr": chr}, {"pos": pos}, {"reference": ref}]},
                        {"alternate": 1, "rsid": 1, "caf": 1})  # projection on desired fields

                    for record in mongo_cursor:
                        record_alts = record["alternate"].split(",")
                        alt_pos = record_alts.index(alt) if alt in record_alts else -1

                        if alt_pos != -1:
                            rsid_list.append(record["rsid"])

                            allelic_freq_list = record["caf"].split(",")
                            if len(record_alts) < len(allelic_freq_list) and ref_freq == "." and alt_freq == ".":
                                ref_freq = allelic_freq_list[0]
                                alt_freq = allelic_freq_list[alt_pos + 1]
                                if len(rsid_list) > 1:
                                    rsid_list.insert(0, rsid_list.pop(-1))

                    if len(rsid_list) > 0:
                        rsid = ",".join(rsid_list)

                    out.write(hline + "\t" + "\t".join([rsid, ref_freq, alt_freq]) + "\n")

        print("Allelic frequencies added successfully!")
        client.close()
        return outfile
    except pymongo.errors.ConnectionFailure as e:
        sys.stderr.write("Could not connect to server: %s\n" % e)
    except pymongo.errors.OperationFailure as o:
        sys.stderr.write("MongoDB interaction failed ({})\n".format(o))
    except FileNotFoundError:
        sys.stderr.write("File {} not found or file error".format(human_file))
    return -1


if __name__ == "__main__":
    import argparse
    import time
    from collections import namedtuple

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ip", default='192.168.10.19', help="IP address of the MongoDB server")
    parser.add_argument("-p", "--port", required=False, help="Port number of the listening MongoDB server",
                        default=27017, type=int)
    parser.add_argument("-v", "--ver", required=True, help="Version of the dbSNP database", )
    parser.add_argument("-u", "--human", required=True, help="Human file to be annotated")
    parser.add_argument("--genome-version", default='hg19', help="Human file to be annotated")

    args = parser.parse_args()
    MongoDB_Info = namedtuple("MongoDB_Info", "ip, port, dbsnp_ver, genome_version")
    mongodb_info = MongoDB_Info(ip=args.ip, port=args.port, dbsnp_ver=args.ver, genome_version=args.genome_version)
    # mongodb_info = MongoDB_Info(ip="localhost", port=27017, dbsnp_ver="150_NEW")

    start = time.perf_counter()
    dbsnp_calculator(args.human, mongodb_info)
    end = time.perf_counter()
    print("--- Elapsed time: {:.2f} seconds ---".format(end - start))

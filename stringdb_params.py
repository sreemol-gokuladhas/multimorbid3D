API_URL = "https://string-db.org/api"
OUTPUT_FORMAT = "tsv"
GET_IDS = "get_string_ids"
INTERACTION_PARTNERS = "interaction_partners"

PARAMS = {
    "species": 9606,
    "limit": 5,
    "echo_query": 1,
    "caller_identity": ""}

def get_params(gene_list):
    params = PARAMS
    params["identifiers"] = "\r".join([str(gene) for gene in gene_list])
    return params

def get_request_ids_url():
    return "/".join([API_URL, OUTPUT_FORMAT, GET_IDS])


def get_request_interaction_url():
    return "/".join([API_URL, OUTPUT_FORMAT, INTERACTION_PARTNERS])

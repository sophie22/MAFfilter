import json
import requests

def build_url(variant: str):
    """ Builds external url path with parameters

    Args:
        variant in the format of Chrom:Pos:Ref:Alt

    Returns:
        url: URL for the API call
    """

    # if len(variant_list) > 1:
    #     variant_list = ",".join(variant_list)
    # else:
    #     variant_list = variant_list[0]

    # ext_url = "%3A".join(variant.split(":"))
    url = "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/{}/annotation?assembly=grch37&include=populationFrequencies&limit=-1&skip=-1&skipCount=false&count=false&Output%20format=json&normalize=false&phased=false&useCache=false&imprecise=false&svExtraPadding=0&cnvExtraPadding=0".format(variant)

    return url

def get_response(url: str):
    """ Make an API query

        Args:
            url (str): Full url to use for the API call. Defaults to None.

        Returns:
            dict: Data from the API call
    """
    for i in range(0, 5):
        try:
            request = requests.get(url, headers={"Accept": "application/json"})
        except Exception as e:
            print("Something went wrong: {}".format(e))
        else:
            if request.ok:
                data = json.loads(request.content.decode("utf-8"))
                if len(data["response"][0]["result"]) > 0:
                    return data["response"][0]["result"][0]
                else:
                #     print("Error {} for URL: {}".format(data["errors"], url))
                    return None
            else:
                print("Error {} for URL: {}".format(request.status_code, url))
                return None

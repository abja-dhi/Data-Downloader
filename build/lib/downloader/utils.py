import os

class utils():
    def __init__(self):
        pass

    def mkch(path):
        try:
            os.mkdir(path)
        except:
            pass
        os.chdir(path)
        return

    def create_log(path, reset, log_name="log.txt"):
        if reset:
            log = open(log_name, "w")
            log.close()
        else:
            if log_name not in os.listdir(path):
                log = open(log_name, "w")
                log.close()
        return
    
    def write_log(error, log_name="log.txt"):
        log = open(log_name, "a")
        if "\n" not in error:
            error = error + "\n"
        log.write(error)
        log.close()

        return
    
class Point:
    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat

    def __str__(self) -> str:
        return str(round(self.lon, 3)) + "_" + str(round(self.lat, 3))
    

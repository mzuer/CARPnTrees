#####################################################################################
########## Script containing definition of classes  used in CARP programm ###########
# Marie Zufferey - UNIL - December 2015 #############################################
# License: Open Source AGPL v3 ######################################################
#####################################################################################

class GenData:
    """ Class that defines an object containing recent polymorphism data (parent class for StrData and SeqData) """

    def __init__(self, data, analysis):
        self.analysis = analysis
        self.data = data
    def __repr__(self):
        return "Object that contains polymorphism data."


class StrData(GenData):
    """ Class that defines an object containing recent polymorphism data in the form of STR data """
    def __init__(self, data):
        GenData.__init__(self, data, "str")
        self.markers = list(data)          # markers name (= columns names of the DF)
        self.labels = list(data.index)     # individuals name (= rows names of the DF)
    def __repr__(self):
        return(super().__repr__()+\
         "\nMarkers:"+str(self.markers)+
        "\nIndiviuals of the data set:"+
        str(self.labels))


class SeqData(GenData):
    """ Class that defines an object containing polymorphism data in the form of sequence data """
    def __init__(self, data):
        import re
        GenData.__init__(self, data, "seq")
        # self.labels = [re.search("(^[A-Z]+[0-9]+)", x.id).group(1) for x in data] # if you want to cut ID name
        for x in range(len(data)):
            # data[x].id = re.search("(^[A-Z]+[0-9]+)", data[x].id).group(1) # if you want to cut ID name
            data[x].id = data[x].id
        self.labels = [x.id for x in data]

    def __repr__(self):
        return(super().__repr__()+\
        "\nIndiviuals of the data set:"+
        str(self.labels))

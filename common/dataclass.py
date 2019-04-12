class DataClass:
    def __init__(self, *variables):
        self.data = []
        for var in variables:
            self.data.append(var)


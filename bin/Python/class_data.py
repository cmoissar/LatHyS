class Data:
    def __init__(self,XYplane=0,XZplane=0,YZplane=0,time=0):
        self.time = time
        self.XY = XYplane
        self.XZ = XZplane
        self.YZ = YZplane
#        self.r=realpart
#        self.i=imagpart

    def vider(self):
        self.time=0
        
class Plane:
    def __init__(self,time=0):
        self.time = time

            

# KernelEnvironment is used for Kernel objects
import numpy as np
import pyopencl as cl

def humanReadableSize(size):
    """
    Converts a size in bytes to human readable form.

    @type size: Number
    @param size: Size in bytes

    @rtype: string
    @return: The size in human readable form.
    """
    if size<=0 or size >=10e18:
        return "%.3g" % size +" B"
    else:
        p=int(np.log(size)/np.log(1024))
        names=["","ki","Mi","Gi","Ti","Pi","Ei"]
        return "%.3g" % (size/1024.0**p)+" "+names[p]+"B"
    
class PIMCKernel:
    def getStats(self):
        """
        Returns information about memory usage of the loaded PIMCKernel.
        Global memory usage is just estimated based on RunParameters and might
        thus not be completely accurate.

        @rtype:   string
        @return:  The information in human readable text format
        """

        ret=""
        dev=self.prg.get_info(cl.program_info.DEVICES)
        if len(dev)!=1:
            raise Exception("Expected the number of devices to be 1, it was "+str(len(dev)))
        dev=dev[0]
        usedGlobalMemory=0
        usedGlobalMemory+=(self.nbrOfThreads+1)*4*4#seeds
        usedGlobalMemory+=(self.nbrOfThreads)*4#accepts
        usedGlobalMemory+=self.runParams.nbrOfWalkers*self.runParams.N*self.system.DOF*4#path
        usedGlobalMemory+=self.nbrOfThreads*self.nbrOfOperators*4#operator
        if self.runParams.enableGlobalOldPath:
            usedGlobalMemory+=self.nbrOfThreads*(2**self.runParams.S-1)*self.system.DOF*4
        if self.runParams.enableBins:
            usedGlobalMemory+=self.runParams.binsPerPart**self.system.DOF*4
            
        ret+=("Global memory (used/max): " +
            humanReadableSize(usedGlobalMemory)+" / "+
            humanReadableSize(dev.get_info(cl.device_info.GLOBAL_MEM_SIZE))+"\n")
        
        ret+=("Local memory (used/max): " +
            humanReadableSize(self.prg.metropolis.get_work_group_info(cl.kernel_work_group_info.LOCAL_MEM_SIZE,dev))+
            " / "+
            humanReadableSize(dev.get_info(cl.device_info.LOCAL_MEM_SIZE))+"\n")
        if self.runParams.enableParallelizePath:
            ret+=("Workgroup size (used/max): "+
                str(self.runParams.N / (2 ** self.runParams.S)*self.runParams.nbrOfWalkersPerWorkGroup)+" / "+
                str(dev.get_info(cl.device_info.MAX_WORK_GROUP_SIZE))+"\n")
        ret+=("Workgroup dimensions (used/max): "+
            str(self.localSize)+" / "+
            str(dev.get_info(cl.device_info.MAX_WORK_ITEM_SIZES))+"\n")
        ret+=("Number of workgroups (used): "+
            str(self.runParams.nbrOfWalkers/self.runParams.nbrOfWalkersPerWorkGroup))    
                
                
        #self.prg.metropolis.get_work_group_info(kernel_work_group_info.LOCAL_MEM_SIZE,dev))
            #CL_KERNEL_WORK_GROUP_SIZE

        return ret

    def getBuildLog(self):
        pass #awaiting mail-list response...Got response, bug fixed, awaiting update...
        """dev=self.prg.get_info(cl.program_info.DEVICES)
        if len(dev)!=1:
            raise Exception("Expected the number of devices to be 1, it was "+str(len(dev)))
        dev=dev[0]
        #return self.prg.get_build_info(dev,cl.program_build_info.LOG)
        return self.prg.get_build_info(dev,cl.program_build_info.LOG)"""

from cmsmltools import *
import os
import tensorflow as tf

out_dir = '/home/users/damante/CMSSW_11_2_1_Patatrack/src/graph_model/'
model_dir = '/home/users/damante/CMSSW_11_2_1_Patatrack/src/model/model_3D4CNN14CNN2_0p00Dropout_0p0010LearningRate/model_3D4CNN14CNN2_0p00Dropout_0p0010LearningRate'

model = tf.keras.models.load_model(model_dir)
print("model successfully loaded")
#final_path = GetDataSetPath(machine)+'graph_model'
#save_graph(final_path, model, variables_to_constants=True)

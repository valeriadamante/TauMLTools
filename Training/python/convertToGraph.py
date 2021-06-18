import cmsmltools
import os
import tensorflow as tf

out_dir = '/home/users/damante/CMSSW_11_2_1_Patatrack/src/graph_model/graph_saved_model.pb'
model_dir = '/home/users/damante/CMSSW_11_2_1_Patatrack/src/model/model_3D4CNN14CNN2_0p00Dropout_0p0010LearningRate/model_3D4CNN14CNN2_0p00Dropout_0p0010LearningRate'


model = tf.keras.models.load_model(model_dir)
print("model successfully loaded")
cmsmltools.save_graph(out_dir, model, variables_to_constants=True)

import torch
import torch.nn as nn
import torchvision.datasets as datasets
import torchvision.transforms as transforms
from torch.autograd import Variable
from torch.utils import data

import scipy.io as scio
import numpy as np



input_size = 11 # low dimension resolution
hidden_size = 100 #
num_classes = 10
num_epochs = 30
batch_size = 10
learning_rate = 0.001

class ShallowNetwork(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        super(ShallowNetwork, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, num_classes)

    def forward(self,x):
        out = self.fc1(x)
        out = self.relu(out)
        out =self.fc2(out)
        return out

# Creating TrainLoader
# load low_dim training data
low_dim_data = scio.loadmat('low_dim_6000.mat')
low_dim = low_dim_data['low_dim'].reshape(6000)
low_dim = np.array([i.reshape(11) for i in low_dim])
low_dim -= low_dim.min()
low_dim /= low_dim.max()
# load labels
orig_data = scio.loadmat('subset_6000.mat')
labels = orig_data['subset_label_X']
labels = labels.reshape(6000)
# convert to dataLoader

tensor_low_dim = torch.Tensor(low_dim[:5000,:])
tensor_labels = torch.Tensor(labels[:5000])
low_dim_dataset = data.TensorDataset(tensor_low_dim, tensor_labels)
train_loader = data.DataLoader(dataset=low_dim_dataset,
                              batch_size=batch_size,)

# Creating TestLoader
# load_low_dim test data
low_dim_data = scio.loadmat('low_dim_test_1000.mat')
low_dim = low_dim_data['low_dim_test'].reshape(1000)
low_dim = np.array([i.reshape(11) for i in low_dim])
low_dim -= low_dim.min()
low_dim /= low_dim.max()
# load labels
labels = low_dim_data['subset_label_X_test']
labels = labels.reshape(1000)
# # convert to dataLoader
tensor_low_dim = torch.Tensor(low_dim[:,:])
tensor_labels = torch.Tensor(labels[:])
low_dim_dataset = data.TensorDataset(tensor_low_dim, tensor_labels)
test_loader = data.DataLoader(dataset=low_dim_dataset,
                              batch_size=batch_size,)

net = ShallowNetwork(input_size, hidden_size, num_classes)
# UNCOMMENT THE NEXT LINE TO USE GPU
# net = net.cuda()
criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(net.parameters(), lr=learning_rate)

for epoch in range(num_epochs):
    for i, (images, labels) in enumerate(train_loader):  
        images = Variable(images)         
        labels = Variable(labels.long())
        
        optimizer.zero_grad()                             # Intialize the hidden weight to all zeros
        outputs = net(images)                           # Forward pass: compute the output class given a image
        loss = criterion(outputs, labels)                 # Compute the loss: difference between the output class and the pre-given label
        loss.backward()                                   # Backward pass: compute the weight
        optimizer.step()                                  # Optimizer: update the weights of hidden nodes
        if (i+1) % 100 == 0:                              # Logging
            print('Epoch [%d/%d], Step [%d/%d], Loss: %.4f'
                 %(epoch+1, num_epochs, i+1, 6000//batch_size, loss.data))

correct = 0
total = 0
for images, labels in test_loader:
    labels = labels.long()
    images = Variable(images)
    outputs = net(images)
    _, predicted = torch.max(outputs.data, 1)  # Choose the best class from the output: The class with the best score
    total += labels.size(0)                    # Increment the total count
    correct += (predicted == labels).sum()     # Increment the correct count
    
print('Accuracy of the network on the 10K test images: %d %%' % (100 * correct / total))
torch.save(net.state_dict(), 'fnn_model.pkl')


# visualize the network
from tensorboardX import SummaryWriter
writer = SummaryWriter('runs/mnist')
writer.add_graph(net,tensor_low_dim)
writer.close()
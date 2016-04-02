from django import forms
from .models import Semester, Mask
from django.conf import settings
import os
import re

# Create the form to hold the mask information
class MaskModelForm(forms.ModelForm):
    class Meta:
        model = Mask
        fields = '__all__'


# Create the form class.
class MaskUploadForm(forms.Form):

    # Create a dropdown based on the run list
    semester = Semester()
    pi_name = forms.ChoiceField(choices=semester.get_allocated_pi(),
                                widget=forms.Select(attrs={'class':'form-control'}))
    file = forms.FileField()

    def parse_maskfile(self, data, pi_index):
        semester = Semester()
        pi_names = [x for (i, x) in semester.get_allocated_pi()]

        # Read the parameters into the data_dict
        data_dict = {}
        for line in data:
            line = line.strip()
            if line != '':
                sline = re.split(b'\t+', line)
                if sline[0] == b'ra':
                    data_dict['ra'] = sline[1].decode('ascii')
                if sline[0] == b'dec':
                    data_dict['dec'] = sline[1].decode('ascii')
                if sline[0] == b'pa':
                    data_dict['position_angle'] = float(sline[1].decode('ascii'))
                if sline[0] == b'label':
                    data_dict['mask'] = sline[1].decode('ascii')
        data_dict['program_id'] = pi_names[pi_index]
        data_dict['trimester'] = semester.get_trimester()

        return data_dict

    def clean_file(self):
        data = self.cleaned_data["file"]
        pi_index = int(self.cleaned_data["pi_name"])
        data_dict = self.parse_maskfile(data, pi_index)

        # Check to see if this mask already exists in the database
        if Mask.objects.filter(mask=data_dict['mask']).exists():
            raise forms.ValidationError(
                "A mask named %s already exists. Mask names *must* be unique." % data_dict['mask'])

        # TODO: We need to save the mask file at some point
        filepath = os.path.join(data_dict['trimester'], data_dict['mask']+'.msk')
        rootname = settings.ROOT_PATH
        upload_path = os.path.join(rootname, 'uploads')
        upload_file = os.path.join(upload_path, filepath)
        f = open(upload_file, 'w')
        for line in data:
            f.write(line.decode('ascii'))
        f.close()


        form = MaskModelForm(data_dict)
        if form.is_valid():
            # Don't put it into the database at this point
            self.instance = form.save(commit=False)
        else:
            # TODO Add more detailed error message
            raise forms.ValidationError(u"The file contains invalid data.")
        return data

    def save(self):
        # We are not overwriting the save method here because form.Form does not have one!
        # We are adding it for convience
        instance = getattr(self, "instance", None)
        if instance:
            instance.save()
        return instance

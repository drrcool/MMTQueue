# -*- coding: utf-8 -*-
# Generated by Django 1.9.5 on 2016-04-06 16:14
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='observing_block',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('duration', models.IntegerField()),
                ('n_visits', models.IntegerField()),
                ('objID', models.CharField(max_length=50)),
                ('start_time', models.CharField(max_length=50)),
                ('ra', models.CharField(max_length=50)),
                ('dec', models.CharField(max_length=50)),
                ('pmra', models.DecimalField(decimal_places=3, max_digits=5)),
                ('pmdec', models.DecimalField(decimal_places=2, max_digits=4)),
                ('pa', models.DecimalField(decimal_places=2, max_digits=6)),
                ('moon', models.CharField(max_length=6)),
                ('exptime', models.DecimalField(decimal_places=2, max_digits=9)),
                ('nexp', models.IntegerField()),
                ('repeats', models.IntegerField()),
                ('oneper', models.BooleanField()),
                ('priority', models.DecimalField(decimal_places=1, max_digits=5)),
                ('grism', models.CharField(max_length=20)),
                ('filters', models.CharField(max_length=20)),
                ('epoch', models.CharField(max_length=20)),
                ('dithersize', models.DecimalField(decimal_places=1, max_digits=5)),
                ('seeing', models.DecimalField(decimal_places=1, max_digits=3)),
                ('photometric', models.BooleanField()),
                ('readtab', models.CharField(max_length=20)),
                ('gain', models.CharField(max_length=20)),
                ('obstype', models.CharField(max_length=20)),
                ('mag', models.DecimalField(decimal_places=2, max_digits=5)),
                ('mask', models.CharField(max_length=20)),
            ],
        ),
    ]
